#include "minimizer.h"
#include <cassert>
#include <cstring>
#include <iostream>
#include "sequence.h"
#include "kvec.h"
#include <numeric>
#include "lchain.h"
/* modified from lh3/minimap2/sketch.c */
/********** start *************/
static inline uint64_t hash64(uint64_t key, uint64_t mask)
{
  key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
  key = key ^ key >> 24;
  key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
  key = key ^ key >> 14;
  key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
  key = key ^ key >> 28;
  key = (key + (key << 31)) & mask;
  return key;
}

typedef struct { // a simplified version of kdq
  int front, count;
  int a[32];
} tiny_queue_t;

static inline void tq_push(tiny_queue_t* q, int x)
{
  q->a[((q->count++) + q->front) & 0x1f] = x;
}

static inline int tq_shift(tiny_queue_t* q)
{
  int x;
  if (q->count == 0) return -1;
  x = q->a[q->front++];
  q->front &= 0x1f;
  --q->count;
  return x;
}

/**
 * Find symmetric (w,k)-minimizers on a DNA sequence
 *
 * @param km     thread-local memory pool; using NULL falls back to malloc()
 * @param str    DNA sequence
 * @param len    length of $str
 * @param w      find a minimizer for every $w consecutive k-mers
 * @param k      k-mer size
 * @param rid    reference ID; will be copied to the output $p array
 * @param is_hpc homopolymer-compressed or not
 * @param p      minimizers
 *               p->a[i].x = kMer<<8 | kmerSpan
 *               p->a[i].y = rid<<32 | lastPos<<1 | strand
 *               where lastPos is the position of the last base of the i-th minimizer,
 *               and strand indicates whether the minimizer comes from the top or the bottom strand.
 *               Callers may want to set "p->n = 0"; otherwise results are appended to p
 */
void mm_sketch(void* km, const char* str, int len, int w, int k, uint32_t rid, int is_hpc, mm128_v* p)
{
  uint64_t shift1 = 2 * (k - 1), mask = (1ULL << 2 * k) - 1, kmer[2] = { 0,0 };
  int i, j, l, buf_pos, min_pos, kmer_span = 0;
  mm128_t buf[256], min = { UINT64_MAX, UINT64_MAX };
  tiny_queue_t tq;

  assert(len > 0 && (w > 0 && w < 256) && (k > 0 && k <= 28)); // 56 bits for k-mer; could use long k-mers, but 28 enough in practice
  memset(buf, 0xff, w * 16);
  memset(&tq, 0, sizeof(tiny_queue_t));
  kv_resize(mm128_t, km, *p, p->n + len / w);

  for (i = l = buf_pos = min_pos = 0; i < len; ++i) {
    int c = nt4_table[(uint8_t)str[i]];
    mm128_t info = { UINT64_MAX, UINT64_MAX };
    if (c < 4) { // not an ambiguous base
      int z;
      if (is_hpc) {
        int skip_len = 1;
        if (i + 1 < len && nt4_table[(uint8_t)str[i + 1]] == c) {
          for (skip_len = 2; i + skip_len < len; ++skip_len)
            if (nt4_table[(uint8_t)str[i + skip_len]] != c)
              break;
          i += skip_len - 1; // put $i at the end of the current homopolymer run
        }
        tq_push(&tq, skip_len);
        kmer_span += skip_len;
        if (tq.count > k) kmer_span -= tq_shift(&tq);
      }
      else kmer_span = l + 1 < k ? l + 1 : k;
      kmer[0] = (kmer[0] << 2 | c) & mask;           // forward k-mer
      // std::cerr << i << " " << kmer[0] << "\n";
      kmer[1] = (kmer[1] >> 2) | (3ULL ^ c) << shift1; // reverse k-mer

      // if (kmer[0] == kmer[1]) continue; // skip "symmetric k-mers" as we don't know it strand
      // z = kmer[0] < kmer[1] ? 0 : 1; // strand
      z = 0; // Do not take strand into account  

      ++l;
      if (l >= k && kmer_span < 256) {
        // std::cerr << rid << " " << z << " " << kmer[z] << " " << i << "\n";
        info.x = hash64(kmer[z], mask) << 8 | kmer_span;
        // std::cerr << info.x << "\n";
        info.y = (uint64_t)rid << 32 | (uint32_t)i << 1 | z;
      }
    }
    else l = 0, tq.count = tq.front = 0, kmer_span = 0;
    buf[buf_pos] = info; // need to do this here as appropriate buf_pos and buf[buf_pos] are needed below
    if (l == w + k - 1 && min.x != UINT64_MAX) { // special case for the first window - because identical k-mers are not stored yet
      for (j = buf_pos + 1; j < w; ++j)
        if (min.x == buf[j].x && buf[j].y != min.y) kv_push(mm128_t, km, *p, buf[j]);
      for (j = 0; j < buf_pos; ++j)
        if (min.x == buf[j].x && buf[j].y != min.y) kv_push(mm128_t, km, *p, buf[j]);
    }
    if (info.x <= min.x) { // a new minimum; then write the old min
      if (l >= w + k && min.x != UINT64_MAX) kv_push(mm128_t, km, *p, min);
      min = info, min_pos = buf_pos;
    }
    else if (buf_pos == min_pos) { // old min has moved outside the window
      if (l >= w + k - 1 && min.x != UINT64_MAX) kv_push(mm128_t, km, *p, min);
      for (j = buf_pos + 1, min.x = UINT64_MAX; j < w; ++j) // the two loops are necessary when there are identical k-mers
        if (min.x >= buf[j].x) min = buf[j], min_pos = j; // >= is important s.t. min is always the closest k-mer
      for (j = 0; j <= buf_pos; ++j)
        if (min.x >= buf[j].x) min = buf[j], min_pos = j;
      if (l >= w + k - 1 && min.x != UINT64_MAX) { // write identical k-mers
        for (j = buf_pos + 1; j < w; ++j) // these two loops make sure the output is sorted
          if (min.x == buf[j].x && min.y != buf[j].y) kv_push(mm128_t, km, *p, buf[j]);
        for (j = 0; j <= buf_pos; ++j)
          if (min.x == buf[j].x && min.y != buf[j].y) kv_push(mm128_t, km, *p, buf[j]);
      }
    }
    if (++buf_pos == w) buf_pos = 0;
  }
  if (min.x != UINT64_MAX)
    kv_push(mm128_t, km, *p, min);
}
/************ end *************/

/* modified from yangao07/abPOA/abpoa_seed.c */
#define sort_key_mm128x(a) ((a).x)
KRADIX_SORT_INIT(mm128x, mm128_t, sort_key_mm128x, 8)
/********** start *************/
int minimizer_t::collect_mm(void* km, const std::vector<seq_t>& seqs, para_t* para) {

  mm_h = (int*)malloc((seqs_size + 2) * sizeof(int));
  mm_h[0] = 0; // record mm i
  for (int i = 0; i < seqs.size(); ++i) { // collect minimizers
    // if (abpt->m > 5) mm_aa_sketch(km, seqs[i], seq_lens[i], abpt->w, abpt->k, i, 0, mm);
    // std::cerr << i << " " << para->m << "\n";
    // std::cerr << i << "\n";
    if (para->m <= 5) mm_sketch(km, seqs[i].seq.c_str(), seqs[i].seq.size(), para->mm_w, para->k, i, 0, &mm_v);
    mm_h[i + 1] = mm_v.n;
  }
  return mm_v.n;
}
minimizer_t::minimizer_t(para_t* para, const std::vector<seq_t>& seqs) {
  init(para, seqs);
};

void minimizer_t::init(para_t* para, const std::vector<seq_t>& seqs) {
  km = nullptr;
  mm_v = { 0, 0, nullptr };
  sorted_mm_v = { 0, 0, nullptr };
  mm_h = nullptr;
  seqs_size = seqs.size();
  ord.resize(seqs_size);
  rid_to_ord.resize(seqs_size);
  std::iota(ord.begin(), ord.end(), 0);
  std::iota(rid_to_ord.begin(), rid_to_ord.end(), 0);
  if (para->progressive_poa || para->enable_seeding) {
    km = km_init();
    collect_mm(km, seqs, para);//mm ->Minimizer
    maxl = 0;
    // std::cerr << "111 " << "\n";
    len.resize(seqs_size);
    for (int i = 0; i < seqs_size; i++) {
      len[i] = seqs[i].seq.size();
      maxl = std::max(maxl, len[i]);
    }
    // std::cerr << "collect_mm" << "\n";

    if (para->enable_seeding) {
      for (int i = 0; i < (int)mm_v.n; ++i) kv_push(mm128_t, km, sorted_mm_v, mm_v.a[i]);
      for (int i = 0; i < seqs_size; i++) {
        radix_sort_mm128x(sorted_mm_v.a + mm_h[i], sorted_mm_v.a + mm_h[i + 1]); // sort mm by k-mer hash values
        // for (int j = mm_h[i]; j < mm_h[i + 1]; j++) {
        //   std::cerr << mm_v.a[j].x << " " << ((mm_v.a[j].y >> 1) & 0x7FFFFFFF) << " ";
        // }
        // std::cerr << "\n";
      }
    }
    // std::cerr << "finish collect_mm" << "\n";
  }

};
void minimizer_t::get_guide_tree(para_t* para) {
  if (para->progressive_poa && seqs_size >= 2) {
    // copy mm1 to mm2
    mm128_v mm_tv = { 0, 0, nullptr };
    for (int i = 0; i < (int)mm_v.n; ++i) kv_push(mm128_t, km, mm_tv, mm_v.a[i]);
    // use mm2 to build guide tree
    if (mm_tv.n == 0) {
      // return ord;
      if(para->verbose) std::cerr << "no minimizer" << "\n";
      kfree(km, mm_tv.a);
      return;
    }

    // if (abpt->verbose >= ABPOA_INFO_VERBOSE) fprintf(stderr, "[%s] Building progressive guide tree ... ", __func__);
    size_t i, _i, j; int rid1, rid2;                                          // mm_hit_n: mimizer hits between each two sequences
    // 0: 0
    // 1: 0 1
    // 2: 0 1 2
    int* mm_hit_n = (int*)calloc((seqs_size * (seqs_size + 1)) >> 1, sizeof(int)); //  ...
    // n: 0 1 ... n-1 n
    // 
    // # total mimizers of i: mm_hit_n[(i*(i+1))/2+i]
    // # total hits for i and j (i>j): mm_hit_n[(i*(i+1)/2)+j]
    // std::cerr << mm_tv.n << "\n";
    radix_sort_mm128x(mm_tv.a, mm_tv.a + mm_tv.n); // sort mm by k-mer hash values
    // std::cerr << "finish sort" << "\n";
    uint64_t last_x = mm_tv.a[0].x;
    int* mm_cnt = (int*)malloc(seqs_size * sizeof(int));
    for (_i = 0, i = 1; i < mm_tv.n; ++i) { // collect mm hits
      if (para->verbose && i % 1000 == 0) std::cerr << "[" << i << "/" << mm_tv.n << "]" << "\n";
      if (mm_tv.a[i].x != last_x) {
        // now [_i, i-1] have the same minimizer k-mer
        memset(mm_cnt, 0, seqs_size * sizeof(int));
        for (j = _i; j < i; ++j) {
          // count mm->a[j]
          rid1 = mm_tv.a[j].y >> 32;
          ++mm_cnt[rid1];
          ++mm_hit_n[((rid1 * (rid1 + 1)) >> 1) + rid1];
        }
        for (rid1 = 0; rid1 < seqs_size - 1; ++rid1) {
          for (rid2 = rid1 + 1; rid2 < seqs_size; ++rid2) {
            mm_hit_n[((rid2 * (rid2 + 1)) >> 1) + rid1] += std::min(mm_cnt[rid1], mm_cnt[rid2]);
          }
        }
        // next minimizer
        last_x = mm_tv.a[i].x, _i = i;
      }
    }
    // now [_i, i-1] have the same minimizer k-mer
    memset(mm_cnt, 0, seqs_size * sizeof(int));
    for (j = _i; j < i; ++j) {
      // count mm->a[j]
      rid1 = mm_tv.a[j].y >> 32;
      ++mm_cnt[rid1];
      ++mm_hit_n[((rid1 * (rid1 + 1)) >> 1) + rid1];
    }
    for (rid1 = 0; rid1 < seqs_size - 1; ++rid1) {
      for (rid2 = rid1 + 1; rid2 < seqs_size; ++rid2) {
        mm_hit_n[((rid2 * (rid2 + 1)) >> 1) + rid1] += std::min(mm_cnt[rid1], mm_cnt[rid2]);
      }
    }
    free(mm_cnt);
    // calculate jaccard similarity between each two sequences
    double* jac_sim = (double*)calloc((seqs_size * (seqs_size)), sizeof(double));
    // 0: 
    // 1: 0 
    // 2: 0 1 
    // std::cerr << "flag" << "\n";
    double max_jac = -1.0, jac; int max_i = -1, max_j = -1;
    for (i = 0; i < (size_t)seqs_size; ++i) {
      for (j = 0; j < (size_t)seqs_size; ++j) {
        if (i == j) continue;
        int shared_n = i > j ? mm_hit_n[((i * (i + 1)) >> 1) + j] : mm_hit_n[((j * (j + 1)) >> 1) + i];
        int tot_n = mm_hit_n[((i * (i + 1)) >> 1) + i] + mm_hit_n[((j * (j + 1)) >> 1) + j] - shared_n;
        if (tot_n == 0) jac = 0;
        else if (tot_n < 0) {
          std::cerr << __func__ << "Bug in progressive tree building. (1)" << "\n";
          exit(EXIT_FAILURE);
        }
        else {
          jac = (0.0 + shared_n) / tot_n;
          jac = jac * ((0.51 * len[i] / maxl) + (0.49 * len[j] / maxl));
          // jac = 1 * jac + 0 * ((0.525 * len[i] / maxl) + (0.475 * len[j] / maxl));
        }
        jac_sim[i * seqs_size + j] = jac; // jac_sim[i][j] = jac_sim[i*(i-1)/2 + j]
        if (jac > max_jac) {
          max_jac = jac; max_i = i, max_j = j;
        }
      }
    }

    // std::cerr << max_jac << " " << max_i << " " << max_j << "\n";
    // build guide tree
    // first pick two with the biggest jac (max_i, max_j)
    max_sim.resize(seqs_size);
    int n_in_map = 2; ord[0] = max_i, ord[1] = max_j;
    rid_to_ord[max_i] = 0;rid_to_ord[max_j] = 1;
    max_sim[max_j] = max_i;
    // then, pick one with biggest jac sum with existing sequence in ord
    while (n_in_map < seqs_size) {
      // std::cerr << n_in_map << " " << seqs_size << "\n";
      max_jac = -1.0, max_j = seqs_size;
      for (rid2 = 0; rid2 < seqs_size; ++rid2) {
        jac = 0.0;
        for (i = 0; i < (size_t)n_in_map; ++i) {
          rid1 = ord[i];
          if (rid1 == rid2) { jac = -1.0; break; }
          // jac = std::max(jac, jac_sim[rid2 * seqs_size + rid1]);
          jac += jac_sim[rid1 * seqs_size + rid2];
        }
        if (jac > max_jac) {
          max_jac = jac;
          max_j = rid2;
        }
      }
      max_i = -1;jac = 0.0;
      for (i = 0; i < (size_t)n_in_map; ++i) {
        rid1 = ord[i];
        if (max_j == rid1) { jac = -1.0; continue; }
        // jac = std::max(jac, jac_sim[rid2 * seqs_size + rid1]);
        if (max_i == -1 || jac < jac_sim[rid1 * seqs_size + max_j]) {
          jac = jac_sim[rid1 * seqs_size + max_j];
          max_i = rid1;
        }
      }
      max_sim[max_j] = max_i;
      // std::cerr << jac << " " << max_sim[max_j] << " " << max_j << "\n";
      if (max_i == seqs_size) {
        std::cerr << __func__ << "Bug in progressive tree building. (2)" << "\n";
        exit(EXIT_FAILURE);
      }
      rid_to_ord[max_j] = n_in_map;
      ord[n_in_map++] = max_j;
    }
    // // calculate jaccard similarity between each two sequences
    // double* jac_sim = (double*)calloc((seqs_size * (seqs_size - 1)) >> 1, sizeof(double));
    // // 0: 
    // // 1: 0 
    // // 2: 0 1 
    // double max_jac = -1.0, jac; int max_i = -1, max_j = -1;
    // for (i = 1; i < (size_t)seqs_size; ++i) {
    //   for (j = 0; j < i; ++j) {
    //     int tot_n = mm_hit_n[((i * (i + 1)) >> 1) + i] + mm_hit_n[((j * (j + 1)) >> 1) + j] - mm_hit_n[((i * (i + 1)) >> 1) + j];
    //     if (tot_n == 0) jac = 0;
    //     else if (tot_n < 0) {
    //       std::cerr << __func__ << "Bug in progressive tree building. (1)" << "\n";
    //       exit(EXIT_FAILURE);
    //     }
    //     else jac = (0.0 + mm_hit_n[((i * (i + 1)) >> 1) + j]) / tot_n;
    //     jac_sim[((i * (i - 1)) >> 1) + j] = jac; // jac_sim[i][j] = jac_sim[i*(i-1)/2 + j]
    //     if (jac > max_jac) {
    //       max_jac = jac; max_i = i, max_j = j;
    //     }
    //   }
    // }

    // // std::cerr << max_jac << " " << max_i << " " << max_j << "\n";
    // // build guide tree
    // // first pick two with the biggest jac (max_i, max_j)
    // int n_in_map = 2; ord[0] = max_j, ord[1] = max_i;

    // // then, pick one with biggest jac sum with existing sequence in ord
    // while (n_in_map < seqs_size) {
    //   // std::cerr << n_in_map << " " << seqs_size << "\n";
    //   max_jac = -1.0, max_i = seqs_size;
    //   for (rid1 = 0; rid1 < seqs_size; ++rid1) {
    //     jac = 0.0;
    //     for (i = 0; i < (size_t)n_in_map; ++i) {
    //       rid2 = ord[i];
    //       if (rid1 == rid2) { jac = -1.0; break; }
    //       else if (rid1 > rid2) jac += jac_sim[((rid1 * (rid1 - 1)) >> 1) + rid2];
    //       else jac += jac_sim[((rid2 * (rid2 - 1)) >> 1) + rid1];
    //     }
    //     if (jac > max_jac) {
    //       max_jac = jac;
    //       max_i = rid1;
    //     }
    //   }
    //   if (max_i == seqs_size) {
    //     std::cerr << __func__ << "Bug in progressive tree building. (2)" << "\n";
    //     exit(EXIT_FAILURE);
    //   }
    //   ord[n_in_map++] = max_i;
    // }

    free(mm_hit_n); free(jac_sim);
    // if (abpt->verbose >= ABPOA_INFO_VERBOSE) fprintf(stderr, "done!\n");

    kfree(km, mm_tv.a);
  }
  // return ord;
};

minimizer_t::~minimizer_t() {
  // std::cerr << "delete" << "\n";
  kfree(km, mm_v.a);kfree(km, sorted_mm_v.a); free(mm_h);
  // std::cerr << 1 << "\n";
  km_destroy(km);
  // std::cerr << 2 << "\n";

}
/*
  minimizers
  mm->a[i].x = kMer<<8 | kmerSpan
  mm->a[i].y = rid<<32 | lastPos<<1 | strand
  anchors
  a[].x: rev<<63 | tid<<32 | tpos
  a[].y: qid << 48 | flags<<40 | q_span<<32 | qpos
*/
int minimizer_t::collect_anchors(mm128_v* anchors, int tid, int qid, int qlen) {
  int i = mm_h[tid], j = mm_h[qid], _i, _j;
  while (i < mm_h[tid + 1] && j < mm_h[qid + 1]) {
    uint64_t xi = sorted_mm_v.a[i].x, xj = sorted_mm_v.a[j].x;
    if (xi == xj) {
      int _i = i;
      for (; _i < mm_h[tid + 1]; ++_i) {
        uint64_t _xi = sorted_mm_v.a[_i].x;
        if (_xi != xi) break;
        uint64_t _yi = sorted_mm_v.a[_i].y;
        for (_j = j; _j < mm_h[qid + 1]; ++_j) {
          uint64_t _xj = sorted_mm_v.a[_j].x;
          if (_xj != xj) break;
          uint64_t _yj = sorted_mm_v.a[_j].y;
          // t_strand<<63 | t_lastPos<<32 | q_lastPos
          mm128_t anchor;
          uint32_t qspan = _xj & 0xff;
          if ((_yi & 1) == (_yj & 1)) { // same strand
            // p->x = (r & 0xffffffff00000000ULL) | rpos;
            // p->y = (uint64_t)q->q_span << 32 | q->q_pos >> 1;
            anchor.x = (_yi & 0x7fffffff00000000ULL) | ((uint32_t)_yi >> 1);
            anchor.y = (uint64_t)qspan << 32 | ((uint32_t)_yj >> 1);
          }
          else { // different strand
            anchor.x = 1ULL << 63 | (_yi & 0xffffffff00000000ULL) | ((uint32_t)_yi >> 1);
            anchor.y = (uint64_t)qspan << 32 | (qlen - (((uint32_t)_yj >> 1) + 1 - qspan) - 1);
          }
          kv_push(mm128_t, km, *anchors, anchor);
        }
      }
      i = _i, j = _j;
    }
    else if (xi < xj) ++i;
    else if (xj < xi) ++j;
  }
  // sort by tpos
  radix_sort_mm128x(anchors->a, anchors->a + anchors->n);
  return anchors->n;
}
/************ end *************/


mm128_t minimizer_t::find_mm(int rid, int tarPos) const {
  mm128_t res = { UINT64_MAX, UINT64_MAX };
  int l = mm_h[rid], r = mm_h[rid + 1];
  while (l + 1 < r) {
    int mid = (l + r) / 2;
    int curPos = (mm_v.a[mid].y >> 1) & 0x7FFFFFFF;
    if (curPos <= tarPos) {
      l = mid;
    }
    else {
      r = mid;
    }
  }
  return mm_v.a[l];
};

mm128_t minimizer_t::find_mm(int& idx, int rid, int tarPos) const {
  mm128_t res = { UINT64_MAX, UINT64_MAX };
  if (tarPos < 0) return res;
  while (idx < mm_h[rid + 1]) {
    int curPos = (mm_v.a[idx].y >> 1) & 0x7FFFFFFF;
    // std::cerr << curPos << " " << mm_v.a[idx].x << "\n";
    if (curPos > tarPos) break;
    if (curPos == tarPos) {
      res = mm_v.a[idx];
    }
    idx++;
  };
  return res;
};
mm128_t minimizer_t::match_mm(uint64_t mm_x, int rid) const {
  // if (mm_x == 2314759) {
  //   for 
  // }
  int l = mm_h[rid], r = mm_h[rid + 1];
  while (l + 1 < r) {
    int mid = (l + r) / 2;
    if (sorted_mm_v.a[mid].x <= mm_x) {
      l = mid;
    }
    else {
      r = mid;
    }
  }
  return sorted_mm_v.a[l];
};

int minimizer_t::dp_chaining(const para_t* para, mm128_v* anchors, int tlen, int qlen) {
  // mg_lchain_dp
  int min_w = para->poa_w + para->k;
  int max_bw = 100, max_dis = 100, max_skip_anchors = 25, max_non_best_anchors = 50, min_local_chain_cnt = 3, min_local_chain_score = 100;
  float chn_pen_gap = 0.8 * 0.01 * para->k, chn_pen_skip = 0.0 * 0.01 * para->k;
  int n_lchains;
  mm128_t* lchains;
  if (para->verbose) std::cerr << "initial anchor size:" << anchors->n << "\n";

  anchors->a = mg_lchain_dp(max_dis, max_dis, max_bw, max_skip_anchors, max_non_best_anchors, min_local_chain_cnt, min_local_chain_score, chn_pen_gap, chn_pen_skip, 0, 1, &anchors->n, anchors->a, &n_lchains, &lchains, km);
  if (para->verbose) std::cerr << "after mg_lchain_dp anchor size:" << anchors->n << "\n";
  if (para->verbose) std::cerr << "after mg_lchain_dp lchains number:" << n_lchains << "\n";
  // sort a by tpos
  radix_sort_mm128x(lchains, lchains + n_lchains);
  // get a complete chain
  // if (para->verbose) {
  //   // std::cerr << (int)(lchains[0].y >> 32) << " " << (int)lchains[0].y << "\n";
  //   for (int i = 0; i < n_lchains; i++) {
  //     uint64_t xi = lchains[i].x, yi = lchains[i].y;
  //     std::cerr << "tpos:" << (xi >> 32) << " " << (int)xi << "\n";
  //     int st = (yi >> 32), ed = (int)yi;
  //     std::cerr << "i [st, ed):"<< i << " " << st << " " << ed << "\n";
  //     for (int i = st; i < ed; i++) {
  //       std::cerr << (int)anchors->a[i].x << " ";
  //     }
  //     std::cerr << "\n";
  //     for (int i = st; i < ed; i++) {
  //       std::cerr << (int)anchors->a[i].y << " ";
  //     }
  //     std::cerr << "\n";
  //   }

  // }
  // find bug

  chain_dp(km, lchains, n_lchains, anchors, min_w, tlen, qlen, max_dis, max_dis, max_bw, chn_pen_gap, chn_pen_skip, 0, 1);
  if (para->verbose) std::cerr << "after chain_dp anchor size:" << anchors->n << "\n";
  kfree(km, lchains);
  return 0;
}

mm128_v minimizer_t::collect_anchors_bycons(const para_t* para, int qid, int qlen, const std::string& cons) {
  // collect cons 's minimizer
  if (para->verbose) std::cerr << "collect cons's minimizer" << "\n";
  int n = seqs_size; // cons 's rid
  mm_v.n = mm_h[n];
  if (para->m <= 5) mm_sketch(km, cons.c_str(), cons.size(), para->mm_w, para->k, n, 0, &mm_v);
  mm_h[n + 1] = mm_v.n;


  sorted_mm_v.n = mm_h[n];
  for (int i = mm_h[n]; i < (int)mm_v.n; ++i) kv_push(mm128_t, km, sorted_mm_v, mm_v.a[i]);
  radix_sort_mm128x(sorted_mm_v.a + mm_h[n], sorted_mm_v.a + mm_h[n + 1]);

  // std::cerr << "mm span\n" << mm_h[qid] << "\n";
  // for (int i = mm_h[qid]; i < mm_h[qid + 1]; i++) {
  //   std::cerr << (sorted_mm_v.a[i].x & 0xff) << " ";
  // }
  // std::cerr << "\n";
  // collect anchor between _seq and cons
  // if (para->verbose) std::cerr << "collect anchor between cons and qseq" << "\n";
  mm128_v anchors = { 0, 0, 0 };
  collect_anchors(&anchors, n, qid, qlen);

  if (para->verbose) std::cerr << "get a optimal chain through dp" << "\n";
  dp_chaining(para, &anchors, cons.size(), qlen);
  // if (para->verbose) {
  //   std::cerr << para->k << " " << anchors.n << "\n";
  //   for (int i = 0; i < anchors.n; i++) {
  //     std::cerr << i << ":" << (int)anchors.a[i].x << " " << (int)(int)anchors.a[i].y << " ";
  //   }
  //   std::cerr << "\n";
  // }
  return anchors;
}