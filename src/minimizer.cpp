#include "minimizer.h"
#include <cassert>
#include <cstring>
#include <iostream>
#include "sequence.h"
#include "kvec.h"
#include "ksort.h"
#include <numeric>

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

  mm_h[0] = 0; // record mm i
  for (int i = 0; i < seqs.size(); ++i) { // collect minimizers
    // if (abpt->m > 5) mm_aa_sketch(km, seqs[i], seq_lens[i], abpt->w, abpt->k, i, 0, mm);
    // std::cerr << i << " " << para->m << "\n";
    // std::cerr << i << "\n";
    if (para->m <= 5) mm_sketch(km, seqs[i].seq.c_str(), seqs[i].seq.size(), para->w, para->k, i, 0, &mm_v);
    mm_h[i + 1] = mm_v.n;
  }
  return mm_v.n;
}
minimizer_t::minimizer_t(para_t* para, const std::vector<seq_t>& seqs) {
  init(para, seqs);
};

void minimizer_t::init(para_t* para, const std::vector<seq_t>& seqs) {
  km = km_init();
  mm_v = { 0, 0, nullptr };
  seqs_size = seqs.size();
  // std::cerr << "collect_mm" << "\n";
  mm_h = (int*)malloc((seqs_size + 1) * sizeof(int));
  collect_mm(km, seqs, para);//mm ->Minimizer
  // std::cerr << "finish collect_mm" << "\n";
};
std::vector<int> minimizer_t::get_guide_tree(para_t* para) {
  std::vector<int> guide_tree(seqs_size);
  std::iota(guide_tree.begin(), guide_tree.end(), 0);
  if (para->progressive_poa && seqs_size > 2) {
    // copy mm1 to mm2
    mm128_v mm_tv = { 0, 0, nullptr };
    for (int i = 0; i < (int)mm_v.n; ++i) kv_push(mm128_t, km, mm_tv, mm_v.a[i]);
    // use mm2 to build guide tree
    if (mm_tv.n == 0) {
      return guide_tree;
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
      // std::cerr << i << " " << mm_tv.n << '\n';
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
    double* jac_sim = (double*)calloc((seqs_size * (seqs_size - 1)) >> 1, sizeof(double));
    // 0: 
    // 1: 0 
    // 2: 0 1 
    double max_jac = -1.0, jac; int max_i = -1, max_j = -1;
    for (i = 1; i < (size_t)seqs_size; ++i) {
      for (j = 0; j < i; ++j) {
        int tot_n = mm_hit_n[((i * (i + 1)) >> 1) + i] + mm_hit_n[((j * (j + 1)) >> 1) + j] - mm_hit_n[((i * (i + 1)) >> 1) + j];
        if (tot_n == 0) jac = 0;
        else if (tot_n < 0) {
          std::cerr << __func__ << "Bug in progressive tree building. (1)" << "\n";
          exit(EXIT_FAILURE);
        }
        else jac = (0.0 + mm_hit_n[((i * (i + 1)) >> 1) + j]) / tot_n;
        jac_sim[((i * (i - 1)) >> 1) + j] = jac; // jac_sim[i][j] = jac_sim[i*(i-1)/2 + j]
        if (jac > max_jac) {
          max_jac = jac; max_i = i, max_j = j;
        }
      }
    }

    // std::cerr << max_jac << " " << max_i << " " << max_j << "\n";
    // build guide tree
    // first pick two with the biggest jac (max_i, max_j)
    int n_in_map = 2; guide_tree[0] = max_j, guide_tree[1] = max_i;

    // then, pick one with biggest jac sum with existing sequence in guide_tree
    while (n_in_map < seqs_size) {
      // std::cerr << n_in_map << " " << seqs_size << "\n";
      max_jac = -1.0, max_i = seqs_size;
      for (rid1 = 0; rid1 < seqs_size; ++rid1) {
        jac = 0.0;
        for (i = 0; i < (size_t)n_in_map; ++i) {
          rid2 = guide_tree[i];
          if (rid1 == rid2) { jac = -1.0; break; }
          else if (rid1 > rid2) jac += jac_sim[((rid1 * (rid1 - 1)) >> 1) + rid2];
          else jac += jac_sim[((rid2 * (rid2 - 1)) >> 1) + rid1];
        }
        if (jac > max_jac) {
          max_jac = jac;
          max_i = rid1;
        }
      }
      if (max_i == seqs_size) {
        std::cerr << __func__ << "Bug in progressive tree building. (2)" << "\n";
        exit(EXIT_FAILURE);
      }
      guide_tree[n_in_map++] = max_i;
    }

    free(mm_hit_n); free(jac_sim);
    // if (abpt->verbose >= ABPOA_INFO_VERBOSE) fprintf(stderr, "done!\n");

    kfree(km, mm_tv.a);
  }
  return guide_tree;
};

minimizer_t::~minimizer_t() {
  kfree(km, mm_v.a); free(mm_h); km_destroy(km);
}
/************ end *************/
