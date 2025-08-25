#include "align.h"
#include <cassert>
#include <algorithm>
#include <immintrin.h>  // AVX2 指令集头文件
#include <cstring>
#include <chrono>


#if __AVX512BW__
constexpr size_t SIMDTotalBytes = 64;
using reg = __m512i;
#elif defined(__AVX2__)
constexpr size_t SIMDTotalBytes = 32;  // 修改处
using reg = __m256i;               // 类型别名用using更规范
#else
constexpr size_t SIMDTotalBytes = 16;
using reg = __m128i;
#endif
constexpr int reg_size = SIMDTotalBytes / sizeof(int); // number of int in reg

constexpr int INF = 0x3f3f3f3f; // 0x3f3f3f3f
constexpr int NEG_INF = 0xc0c0c0c0; //
constexpr int M_OP = 1, D_OP = 2, I_OP = 4, ALL_OP = 7;

inline int calj(int acj, int Bs) {
  return  acj - Bs * reg_size;
}
std::vector<res_t> POA(para_t* para, const graph& DAG, const std::string& seq) {
  int n = DAG.node.size(), m = seq.size();
  assert(seq[0] == para->m - 1);
  const std::vector<node_t>& node = DAG.node;const std::vector<int>& rank = DAG.rank;
  std::vector<int> mat = para->mat;
  int match = para->match, mismatch = para->mismatch, para_m = para->m, e1 = para->gap_ext1, o1 = para->gap_open1 + e1;
  // std::cout << mat << " " << mis << " " << o1 << " " << e1 << "\n";
  std::vector<std::vector<int>> M(n, std::vector<int>(m + 1, NEG_INF)), D(M), I(M);
  M[0][0] = 0;
  // std::cout << "n:" << n << " " << "m:" << m << "\n";

  for (int i = 1; i < n; i++) {
    const node_t& cur = node[rank[i]];
    // if (i % 100 == 0) std::cout << i << "\n";
    for (int k = 0; k < cur.in.size(); k++) {
      int p = node[cur.in[k]].rank; // rank
      const node_t& pre = node[cur.in[k]];
      for (int j = 1; j <= m; j++) {
        if (j - 1 >= 0) M[i][j] = std::max(M[i][j], M[p][j - 1] + (pre.base == seq[j - 1] ? match : mismatch)); // qsource mat[pre.base * para_m + seq[j - 1]]
        D[i][j] = std::max({ D[i][j], D[p][j] + e1, M[p][j] + o1 });    // dsource
      }
    }
    for (int j = 1; j <= m; j++) {
      if (j - 1 >= 0) I[i][j] = std::max(I[i][j - 1] + e1, M[i][j - 1] + o1); // isource
      M[i][j] = std::max({ M[i][j], D[i][j], I[i][j] });  // three source
    }
  }
  // std::cout << M[n - 1][m] << "\n";
  // std::cout << "M:\n";
  // for (int i = 0; i < n; i++) {
  //   for (int j = 0; j <= m; j++) {
  //     std::cout << M[i][j] << " \n"[j == m];
  //   }
  // }
  // std::cout << "D:\n";
  // for (int i = 0; i < n; i++) {
  //   for (int j = 0; j <= m; j++) {
  //     std::cout << D[i][j] << " \n"[j == m];
  //   }
  // }
  // std::cout << "I:\n";
  // for (int i = 0; i < n; i++) {
  //   for (int j = 0; j <= m; j++) {
  //     std::cout << I[i][j] << " \n"[j == m];
  //   }
  // }
  int i = n - 1, j = m;
  std::vector<res_t> res;
  char op = 'M';
  // std::cerr << "finsh align" << "\n";
  while (i > 0 || j > 0) {
    // dsource
    if (op == 'M') {
      if (M[i][j] == D[i][j]) op = 'D';
      if (M[i][j] == I[i][j]) op = 'I';
    }
    // std::cerr << op << " " << i << " " << j << "\n";
    if (op == 'M') {
      const node_t& cur = node[rank[i]];
      for (int k = 0; k < cur.in.size(); k++) {
        int p = node[cur.in[k]].rank; // rank
        const node_t& pre = node[cur.in[k]];
        if (j - 1 >= 0 && M[i][j] == M[p][j - 1] + (pre.base == seq[j - 1] ? match : mismatch)) {
          if (pre.base == seq[j - 1]) {
            // std::cout << "M";
            res.emplace_back(res_t(pre.id, pre.base));
          }
          else {
            // std::cout << "X";
            const node_t& par = node[pre.par_id]; // dsu.find par
            if (par.aligned_node[seq[j - 1]] != -1) {
              res.emplace_back(res_t(par.aligned_node[seq[j - 1]], seq[j - 1]));
            }
            else res.emplace_back(res_t(-1, seq[j - 1], pre.par_id));
          }
          i = p, j--;
          break;
        }
      }
      continue;
    }
    if (op == 'I') {
      if (j - 1 >= 0) {
        // std::cout << "I";
        res.emplace_back(res_t(-1, seq[j - 1]));
        if (I[i][j] == I[i][j - 1] + e1) {
          j--;
          continue;
        }
        if (I[i][j] == M[i][j - 1] + o1) {
          op = 'M';
          j--;
          continue;
        }
      }
    }
    if (op == 'D') {
      const node_t& cur = node[rank[i]];
      // std::cerr << M[i][j] << " " << D[i][j] << "\n";
      for (int k = 0; k < cur.in.size(); k++) {
        int p = node[cur.in[k]].rank; // rank
        const node_t& pre = node[cur.in[k]];
        if (D[i][j] == D[p][j] + e1) {
          i = p;
          break;
        }
        if (D[i][j] == M[p][j] + o1) {
          op = 'M';
          i = p;
          break;
        }
      }
      continue;
    }
  }
  // std::cout << "sorce:" << M[n - 1][m] << "\n";
  return res;
};

std::vector<res_t> abPOA(const para_t* para, const graph* DAG, const std::string& _seq, aligned_buff_t* mpool) {
  // std::chrono::microseconds total_part1(0);
  // std::chrono::microseconds total_part2(0);
  // std::chrono::microseconds total_part3(0);
  // auto start1 = std::chrono::high_resolution_clock::now();

  int n = DAG->node.size(), m = _seq.size();
  std::string seq = _seq;
  assert(seq[0] == para->m - 1);
  const std::vector<node_t>& node = DAG->node;const std::vector<int>& rank = DAG->rank;
  std::vector<int> mat = para->mat;
  int match = para->match, mismatch = para->mismatch, para_m = para->m, e1 = para->gap_ext1, o1 = para->gap_open1 + e1;
  reg MATCH = _mm256_set1_epi32(match), MISMATCH = _mm256_set1_epi32(mismatch), E1 = _mm256_set1_epi32(e1), O1 = _mm256_set1_epi32(o1), Neg_inf = _mm256_set1_epi32(NEG_INF);
  // std::cout << mat << " " << mis << " " << o1 << " " << e1 << "\n";
  int col_size = (m + 1 + reg_size - 1) / reg_size * reg_size;
  seq.append(col_size - m, char26_table['N']);
  // std::cerr << seq.size() << " " << col_size << "\n";
  assert(col_size % reg_size == 0);
  std::vector<int> Ms(n, m + 1), Me(n); // index is rank  [Ms[i], Me[i]]
  std::vector<int> Bs(n), Be(n);  // index is rank    [Bs[i], Be[i] - 1)  Be[i] - 1 's block is Neg_inf is designed for M direciton dp
  size_t mtx_size = 0;
  const int w = para->b + m / para->f;
  size_t sum = 0;
  for (int i = 0; i < n; i++) {
    Ms[i] = std::max(0, std::min(DAG->hlen[i], m - DAG->tlen[i]) - w), Me[i] = std::min(m, std::max(DAG->hlen[i], m - DAG->tlen[i]) + w);
    Bs[i] = Ms[i] / reg_size, Be[i] = Me[i] / reg_size + 2; // [block_s,block_e)
    mtx_size += (Be[i] - Bs[i]) * reg_size;
    // sum += Me[i] - Ms[i] + 1;
  }
  void* buff = nullptr;
  // ===== 2. 对齐内存分配 =====
  // if ((ret = posix_memalign(&buff, SIMDTotalBytes, 3 * mtx_size * sizeof(int))) != 0) {
  // std::cerr << "ac_mtx_size:" << 1ll * n * m << "\n";
  // std::cerr << "mtx_size:" << sum << "\n";
  // std::cerr << m << " " << sum * 1.0 / n << "\n";
  // std::cerr << sum * 1.0 / (1ll * (n * m)) << "\n";
  // return std::vector<res_t>();
  // std::cerr << buff << "\n";
  if (mpool != nullptr) mpool->alloc_aligned(&buff, SIMDTotalBytes, 3 * mtx_size * sizeof(int));
  else alloc_aligned(&buff, SIMDTotalBytes, 3 * mtx_size * sizeof(int)); // malloc
  std::vector<int*> M(n), D(n), I(n);
  M[0] = (int*)buff;
  D[0] = M[0] + mtx_size;
  I[0] = D[0] + mtx_size;
  for (int i = 1; i < n; i++) {
    int offset = (Be[i - 1] - Bs[i - 1]) * reg_size;
    M[i] = M[i - 1] + offset;
    D[i] = D[i - 1] + offset;
    I[i] = I[i - 1] + offset;
  }
  // std::cerr << "n:" << n << " " << w << "\n";
  // std::vector<int> R = DAG->calculateR(); // index is rank
  // std::cerr << "R:" << node[1].rank << " " << R[29650] << "\n";
  // Ms[0] = Me[0] = 0;

  // Bs[0] = std::max(0, std::min(Ms[0], m - R[0]) - w), Be[0] = std::min(m, std::max(Me[0], m - R[0]) + w); // [Bs,Be]
  // int tbs = std::max(0, std::min(DAG->hlen[0], m - DAG->tlen[0]) - w), tbe = std::min(m, std::max(DAG->hlen[0], m - DAG->tlen[0]) + w);
  // Bs[0] = tbs, Be[0] = tbe;

  // return std::vector<res_t>();
  // int block_s = Bs[0] / reg_size, block_e = Be[0] / reg_size + 1; // [block_s,block_e)
  for (int bid = Bs[0]; bid < Be[0]; bid++) {
    // memset(M, 0xc0, col_size * sizeof(int));//M[0] = NEG_INF
    // memset(D, 0xc0, col_size * sizeof(int));//D[0] = NEG_INF
    int* M_i = M[0];
    int* D_i = D[0];
    _mm256_store_si256((reg*)(M_i + bid * reg_size), Neg_inf);
    _mm256_store_si256((reg*)(D_i + bid * reg_size), Neg_inf);
  }
  M[0][0] = 0;
  // if (block_s - 1 >= 0) _mm256_store_si256((reg*)(M + (block_s - 1) * reg_size), Neg_inf);
  // if (block_e < col_size / reg_size) _mm256_store_si256((reg*)(M + block_e * reg_size), Neg_inf);
  // memset(M, 0xc0, col_size * sizeof(int));//M[0] = NEG_INF
  // memset(D, 0xc0, col_size * sizeof(int));//D[0] = NEG_INF
  // memset(buff, 0xc0, 3 * mtx_size * sizeof(int)); // 初始化为-INF 0xc0c0c0c0
  // for (int k = 0; k < node[0].out.size(); k++) {
  //   int suc = node[node[0].out[k]].rank;
  //   Ms[suc] = std::min(Ms[suc], 0 + 1), Me[suc] = std::max(Me[suc], 0 + 1);
  // }
  // auto end1 = std::chrono::high_resolution_clock::now();
  // total_part1 += std::chrono::duration_cast<std::chrono::microseconds>(end1 - start1);
  int tot = 0;
  // std::cerr << "node:" << node.size() << "\n";
  // std::cerr << "path len:" << DAG->hlen[node[1].rank] << " " << DAG->tlen[node[0].rank] << "\n";
  // std::cerr << "m:" << m << "\n";
  for (int i = 1; i < n; i++) {
    const node_t& cur = node[rank[i]];
    int* M_i = M[i];
    int* D_i = D[i];
    int block_num = Be[i] - Bs[i] - 1; // not contain Be[i] - 1 's block
    // tot += cur.in.size();
    // Bs[i] = std::max(0, std::min(Ms[i], m - R[i]) - w), Be[i] = std::min(m, std::max(Me[i], m - R[i]) + w); // [Bs,Be]
    // int tbs = std::max(0, std::min(DAG->hlen[i], m - DAG->tlen[i]) - w), tbe = std::min(m, std::max(DAG->hlen[i], m - DAG->tlen[i]) + w);
    // Bs[i] = tbs, Be[i] = tbe;

    // int tbs = std::max(0, std::min(DAG->hmin[i], m - DAG->tmax[i]) - w), tbe = std::min(m, std::max(DAG->hmax[i], m - DAG->tmin[i]) + w);
    // std::cerr << m << " " << Bs[i] << " " << Be[i] << "\n";
    // block_s = Bs[i] / reg_size, block_e = Be[i] / reg_size + 1; // [block_s,block_e)
    sum += Me[i] - Ms[i] + 1;
    // std::cerr << "B:" << m << " " << w << " " << tbs << " " << tbe << '\n';
    auto start2 = std::chrono::high_resolution_clock::now();
    for (int k = 0; k < cur.in.size(); k++) {
      int p = node[cur.in[k]].rank; // rank
      const node_t& pre = node[cur.in[k]];
      reg PRE_BASE = _mm256_set1_epi32(pre.base);
      int prev = NEG_INF;
      char prech = char26_table['N'];
      int* M_p = M[p];
      int* D_p = D[p];
      for (int bid = 0; bid < block_num; bid++) { // SIMD
        int j = bid * reg_size;
        int acBid = (Bs[i] + bid);
        int acj = acBid * reg_size;
        reg Mij = k == 0 ? Neg_inf : _mm256_load_si256((reg*)(M_i + j));
        reg Dij = k == 0 ? Neg_inf : _mm256_load_si256((reg*)(D_i + j));
        reg SEQ; // seq[j - 1]
        int pj = calj(acj, Bs[p]);
        // std::cerr << "debug" << pj << "\n";
        if (!((std::max(0, acBid - 1) < Bs[pre.rank]) || (acBid > Be[pre.rank] - 2))) { //Bs[pre.rank] <= j - 1 && j <= Be[pre.rank]
          // std::cerr << "debug: " << bid << " " << Bs[pre.rank] / reg_size << " " << Be[pre.rank] / reg_size << "\n";
          if (acj == 0) { // 特殊处理 j=0 的情况
            alignas(32) int tmp[8] = { prech,seq[0], seq[1], seq[2],seq[3], seq[4], seq[5], seq[6] };
            SEQ = _mm256_load_si256((reg*)tmp);
          }
          else {  // 直接加载连续内存 (高效)
            const __m128i seq_chunk = _mm_loadu_si128((const __m128i*)(seq.c_str() + acj - 1));
            SEQ = _mm256_cvtepi8_epi32(seq_chunk);  // 符号扩展到32位
          }
          // M[i][j] = std::max(M[i][j], M[p][j - 1] + mat[pre.base * para_m + seq[j - 1]]); // qsource
          reg TMP = prev == NEG_INF ? _mm256_setr_epi32(prev, M_p[pj], M_p[pj + 1], M_p[pj + 2], M_p[pj + 3], M_p[pj + 4], M_p[pj + 5], M_p[pj + 6]) : _mm256_loadu_si256((reg*)(M_p + pj - 1)); //M[p][j - 1]

          reg mask = _mm256_cmpeq_epi32(PRE_BASE, SEQ);
          // pre.base == seq[j - 1] ? match : mismatch
          TMP = _mm256_add_epi32(TMP, _mm256_blendv_epi8(MISMATCH, MATCH, mask));
          Mij = _mm256_max_epi32(Mij, TMP); //std::max(M[i][j], M[p][j - 1] + mat[pre.base * para_m + seq[j - 1]]); 
        }
        if (Bs[pre.rank] <= acBid && acBid < Be[pre.rank] - 1) { // Bs[pre.rank] <= j && j <= Be[pre.rank]
          // D[i][j] = std::max({ D[i][j], D[p][j] + e1, M[p][j] + o1 });    // dsource
          reg TMP = _mm256_load_si256((reg*)(D_p + pj));
          Dij = _mm256_max_epi32(Dij, _mm256_add_epi32(TMP, E1)); // max(D[i][j], D[p][j] + e1)
          TMP = _mm256_load_si256((reg*)(M_p + pj));
          Dij = _mm256_max_epi32(Dij, _mm256_add_epi32(TMP, O1)); // max(D[i][j], M[p][j] + o1)
        }
        _mm256_store_si256((reg*)(M_i + j), Mij); // M[i][j] = max
        _mm256_store_si256((reg*)(D_i + j), Dij); // D[i][j] = max
        prev = M_p[pj + reg_size - 1]; // pre = M[p][j +reg_size - 1]
        // prech = seq[j + reg_size - 1];// prech =seq[j + reg_size - 1]

      }
    }
    // Be[i] - 1
    // if (block_s - 1 >= 0) _mm256_store_si256((reg*)(M_i + (block_s - 1) * reg_size), Neg_inf);
    _mm256_store_si256((reg*)(M_i + block_num * reg_size), Neg_inf);

    // if (i % 100 == 0) std::cout << i << "\n";
    // auto end2 = std::chrono::high_resolution_clock::now();
    // total_part2 += std::chrono::duration_cast<std::chrono::microseconds>(end2 - start2);
    // auto start3 = std::chrono::high_resolution_clock::now();
    // I[i][0] = NEG_INF;
    int* I_i = I[i];
    // I_i[std::max(1, Bs[i] * reg_size) - 1] = M_i[std::max(1, Bs[i] * reg_size) - 1] = NEG_INF;
    // int max_pos = -1;
    I_i[0] = NEG_INF;
    M_i[0] = std::max({ M_i[0], D_i[0], I_i[0] });
    for (int j = 1; j < block_num * reg_size; j++) { // block_num * reg_size
      I_i[j] = std::max(I_i[j - 1] + e1, M_i[j - 1] + o1); // isource
      M_i[j] = std::max({ M_i[j], D_i[j], I_i[j] });  // three source
      // max_pos = (max_pos == -1 || M_i[j] > M_i[max_pos]) ? j : max_pos;
    }
    // if (debug) {
    //   // for (int j = std::max(1, block_s * reg_size); j < block_e * reg_size; j++) {
    //   //   std::cerr << i << " " << j << " " << M_i[j] << "\n";
    //   // }
    //   return std::vector<res_t>();
    // }
    // auto end3 = std::chrono::high_resolution_clock::now();
    // total_part3 += std::chrono::duration_cast<std::chrono::microseconds>(end3 - start3);
    // for (int k = 0; k < cur.out.size(); k++) {
    //   int suc = node[cur.out[k]].rank;
    //   Ms[suc] = std::min(Ms[suc], max_pos + 1), Me[suc] = std::max(Me[suc], max_pos + 1);
    // }
  }

  // std::cerr << "avg pre size:" << tot / n << "\n";
  // std::cerr << "部分1总耗时: " << total_part1.count() / 1000 << " ms\n";
  // std::cerr << "部分2总耗时: " << total_part2.count() / 1000 << " ms\n";
  // std::cerr << "部分3总耗时: " << total_part3.count() / 1000 << " ms\n";
  int i = n - 1, acj = m;
  int ans = M[i][calj(acj, Bs[i])];
  // std::cout << "M:\n";
  // for (int i = 0; i < n; i++) {
  //   std::cout << "Be:" << Bs[i] << " " << Be[i] << "\n";
  //   for (int j = 0; j < Be[i] * reg_size; j++) {
  //     std::cout << M[i][j] << " ";
  //   }
  //   std::cout << "\n";
  // }
  // std::cout << "D:\n";
  // for (int i = 0; i < n; i++) {
  //   std::cout << "Be:" << Bs[i] << " " << Be[i] << "\n";
  //   for (int j = 0; j < Be[i] * reg_size; j++) {
  //     std::cout << D[i][j] << " ";
  //   }
  //   std::cout << "\n";

  // }
  // std::cout << "I:\n";
  // for (int i = 0; i < n; i++) {
  //   std::cout << "Be:" << Bs[i] << " " << Be[i] << "\n";
  //   for (int j = 0; j < Be[i] * reg_size; j++) {
  //     std::cout << I[i][j] << " ";
  //   }
  //   std::cout << "\n";
  // }
  // return std::vector<res_t>();
  // if (3 * mtx_size == 5208) {
  //   std::cout << "M:\n";
  //   for (int i = 0; i < n; i++) {
  //     std::cout << "Be:" << Bs[i] << " " << Be[i] << "\n";
  //     for (int j = 0; j < col_size; j++) {
  //       std::cout << M[i * col_size + j] << " \n"[j + 1 == col_size];
  //     }
  //   }
  //   std::cout << "D:\n";
  //   for (int i = 0; i < n; i++) {
  //     std::cout << "Be:" << Bs[i] << " " << Be[i] << "\n";
  //     for (int j = 0; j < col_size; j++) {
  //       std::cout << D[i * col_size + j] << " \n"[j + 1 == col_size];
  //     }
  //   }
  //   std::cout << "I:\n";
  //   for (int i = 0; i < n; i++) {
  //     std::cout << "Be:" << Bs[i] << " " << Be[i] << "\n";
  //     for (int j = 0; j < col_size; j++) {
  //       std::cout << I[i * col_size + j] << " \n"[j + 1 == col_size];
  //     }
  //   }
  // }
  // std::cerr << "finish" << "\n";
  // if (3 * mtx_size == 5208) std::cerr << _seq.size() << " " << seq.size() << "\n";
  std::vector<res_t> res;
  // return res;
  int j = calj(acj, Bs[i]);
  int op = D[i][j] >= M[i][j] ? D_OP : I[i][j] >= M[i][j] ? I_OP : M_OP;
  // std::cerr << i << " " << j << " " << M[i][calj(acj, Bs[i])] << "\n";
  // std::cerr << "finsh align" << "\n";
  while (i > 0 || acj > 0) {
    // dsource
    j = calj(acj, Bs[i]);
    // if (ans == 5114 && op == 'M') {
      // std::cerr << M[i][j] << "\n";
      // std::cerr << D[i][j] << "\n";
      // std::cerr << I[i][j] << "\n";
    // }
    // if (ans == 5114) std::cerr << op << " " << i << " " << acj << "\n";
    // if (ans == 5114 && op == 'M' && i == 6840 && acj == 1182) {
    //   // exit(0);
    // }
    // if (3 * mtx_size == 3237179904) {
    //   if (op == 'M') {
    //     std::cerr << M[i * col_size + j] << "\n";
    //   }
    //   else if (op == 'D')
    //   {
    //     std::cerr << D[i * col_size + j] << "\n";
    //   }
    //   if (i == 0 && j == 51) {
    //     return res;
    //   }
    //   if (i == 0) exit(1);
    // }
    // std::cerr << op << "\n";
    // if (op & M_OP) {  // M
    //   // std::cerr << op << " " << i << " " << j << " ";
    //   // std::cerr << M[i][j] << " " << D[i][j] << " " << I[i][j] << "\n";
    //   // if (i == 15932 && j == 205) exit(0);
    //   const node_t& cur = node[rank[i]];
    //   int bk = -1;
    //   for (int k = 0; k < cur.in.size(); k++) {
    //     int p = node[cur.in[k]].rank; // rank
    //     const node_t& pre = node[cur.in[k]];
    //     if (pre.base != seq[acj - 1]) continue;
    //     // M
    //     int pj = calj(acj, Bs[p]);
    //     // if (ans == 5114 && pj - 1 >= 0) 
    //       // std::cerr << M[p][pj - 1] << " " << (pre.base == seq[acj - 1] ? match : mismatch) << "\n";
    //     if (pj - 1 >= 0 && M[i][j] == M[p][pj - 1] + match && (bk == -1 || cur.in_weight[k] > cur.in_weight[bk])) {
    //       bk = k;
    //       // break;
    //     }
    //   }
    //   // std::cerr << bk << "\n";
    //   // assert(bk != -1);
    //   if (bk != -1) {
    //     // std::cerr << "M";
    //     op = ALL_OP;
    //     int p = node[cur.in[bk]].rank; // rank
    //     const node_t& pre = node[cur.in[bk]];
    //     if (pre.base == seq[acj - 1]) {
    //       // std::cout << "M";
    //       res.emplace_back(res_t(pre.id, pre.base));
    //     }
    //     else {
    //       // std::cout << "X";
    //       const node_t& par = node[pre.par_id]; // dsu.find par
    //       if (par.aligned_node[seq[acj - 1]] != -1) {
    //         res.emplace_back(res_t(par.aligned_node[seq[acj - 1]], seq[acj - 1]));
    //       }
    //       else res.emplace_back(res_t(-1, seq[acj - 1], pre.par_id));
    //     }
    //     i = p, acj--;
    //     continue;
    //   }
    // }
    if (op & I_OP) {
      const node_t& cur = node[rank[i]];
      if (j - 1 >= 0) {
        if (!(op & M_OP) || M[i][j] == I[i][j]) {
          // std::cerr << "I";
          res.emplace_back(res_t(-1, seq[acj - 1]));
          if (I[i][j] == I[i][j - 1] + e1) {
            op = I_OP;
            acj--;
            continue;
          }
          if (I[i][j] == M[i][j - 1] + o1) {
            op = M_OP | D_OP;
            acj--;
            continue;
          }
        }
      }
    }
    if (op & D_OP) {
      if (!(op & M_OP) || M[i][j] == D[i][j]) {
        // std::cerr << "D";
        const node_t& cur = node[rank[i]];
        // std::cerr << M[i * col_size + j] << " " << D[i * col_size + j] << "\n";
        int bk = -1;char bop;
        for (int k = 0; k < cur.in.size(); k++) {
          int p = node[cur.in[k]].rank; // rank
          const node_t& pre = node[cur.in[k]];
          if (Bs[pre.rank] <= acj / reg_size && acj / reg_size < Be[pre.rank] - 1) {
            int pj = calj(acj, Bs[p]);
            if (D[i][j] == D[p][pj] + e1 && (bk == -1 || cur.in_weight[k] > cur.in_weight[bk])) {
              bk = k;
              bop = D_OP;
              // break;
            }
            // if (D[i][j] == M[p][pj] + o1 && (bk == -1 || cur.in_weight[k] > cur.in_weight[bk])) {
            //   bk = k;
            //   bop = 'M';
            //   // break;
            // }
          }
        }
        if (bk == -1) {
          for (int k = 0; k < cur.in.size(); k++) {
            int p = node[cur.in[k]].rank; // rank
            const node_t& pre = node[cur.in[k]];
            if (Bs[pre.rank] <= acj / reg_size && acj / reg_size < Be[pre.rank] - 1) {
              int pj = calj(acj, Bs[p]);
              if (D[i][j] == M[p][pj] + o1 && (bk == -1 || cur.in_weight[k] > cur.in_weight[bk])) {
                bk = k;
                bop = M_OP | I_OP;
                // break;
              }
            }
          }
        }
        assert(bk != -1);
        i = node[cur.in[bk]].rank;
        op = bop;
        continue;

      }
    }
    // if (op &) {
    //   if (M[i][j] == D[i][j]) op = 'D';
    //   if (op == 'M' && M[i][j] == I[i][j]) op = 'I';
    // }
    if (op & M_OP) {  // M
      // std::cerr << op << " " << i << " " << j << " ";
      // std::cerr << M[i][j] << " " << D[i][j] << " " << I[i][j] << "\n";
      // if (i == 15932 && j == 205) exit(0);
      const node_t& cur = node[rank[i]];
      int bk = -1;
      for (int k = 0; k < cur.in.size(); k++) {
        int p = node[cur.in[k]].rank; // rank
        const node_t& pre = node[cur.in[k]];
        if (pre.base != seq[acj - 1]) continue;
        // M
        int pj = calj(acj, Bs[p]);
        // if (ans == 5114 && pj - 1 >= 0) 
          // std::cerr << M[p][pj - 1] << " " << (pre.base == seq[acj - 1] ? match : mismatch) << "\n";
        if (pj - 1 >= 0 && M[i][j] == M[p][pj - 1] + match && (bk == -1 || cur.in_weight[k] > cur.in_weight[bk])) {
          bk = k;
          // break;
        }
      }
      // std::cerr << bk << "\n";
      // assert(bk != -1);
      if (bk != -1) {
        // std::cerr << "M";
        op = ALL_OP;
        int p = node[cur.in[bk]].rank; // rank
        const node_t& pre = node[cur.in[bk]];
        if (pre.base == seq[acj - 1]) {
          // std::cout << "M";
          res.emplace_back(res_t(pre.id, pre.base));
        }
        else {
          // std::cout << "X";
          const node_t& par = node[pre.par_id]; // dsu.find par
          if (par.aligned_node[seq[acj - 1]] != -1) {
            res.emplace_back(res_t(par.aligned_node[seq[acj - 1]], seq[acj - 1]));
          }
          else res.emplace_back(res_t(-1, seq[acj - 1], pre.par_id));
        }
        i = p, acj--;
        continue;
      }
    }
    if (op & M_OP) {  // X
      // std::cerr << op << " " << i << " " << j << " ";
      // std::cerr << M[i][j] << " " << D[i][j] << " " << I[i][j] << "\n";
      // if (i == 15932 && j == 205) exit(0);
      const node_t& cur = node[rank[i]];
      int bk = -1;

      for (int k = 0; k < cur.in.size(); k++) {
        int p = node[cur.in[k]].rank; // rank
        const node_t& pre = node[cur.in[k]];
        if (pre.base == seq[acj - 1]) continue;
        int pj = calj(acj, Bs[p]);
        // X
        // if (ans == 5114 && pj - 1 >= 0) 
          // std::cerr << M[p][pj - 1] << " " << (pre.base == seq[acj - 1] ? match : mismatch) << "\n";
        if (pj - 1 >= 0 && M[i][j] == M[p][pj - 1] + mismatch && (bk == -1 || cur.in_weight[k] > cur.in_weight[bk])) {
          bk = k;
          // break;
        }
      }
      // std::cerr << bk << "\n";
      assert(bk != -1);
      if (bk != -1) {
        // std::cerr << "X";
        op = ALL_OP;
        int p = node[cur.in[bk]].rank; // rank
        const node_t& pre = node[cur.in[bk]];
        if (pre.base == seq[acj - 1]) {
          // std::cout << "M";
          res.emplace_back(res_t(pre.id, pre.base));
        }
        else {
          // std::cout << "X";
          const node_t& par = node[pre.par_id]; // dsu.find par
          if (par.aligned_node[seq[acj - 1]] != -1) {
            res.emplace_back(res_t(par.aligned_node[seq[acj - 1]], seq[acj - 1]));
          }
          else res.emplace_back(res_t(-1, seq[acj - 1], pre.par_id));
        }
        i = p, acj--;
        continue;
      }
    }
    // if (op == 'M') {
    //   if (M[i][j] == I[i][j]) op = 'I';
    //   if (op == 'M' && M[i][j] == D[i][j]) op = 'D';
    // }

  }
  // std::cerr << "\n";
  // if (acj > 0) {
  //   std::cerr << "acj:" << acj << "\n";
  //   exit(0);
  // }
  // std::cerr << "finish" << "\n";
  // std::cout << "sorce:" << M[n - 1][m] << "\n";
  if (mpool == nullptr) free_aligned(buff);
  return res;
}

std::vector<res_t> POA_SIMD_ORIGIN(const para_t* para, const graph* DAG, const std::string& _seq, aligned_buff_t* mpool) {
  std::chrono::microseconds total_part1(0);
  std::chrono::microseconds total_part2(0);
  std::chrono::microseconds total_part3(0);
  auto start1 = std::chrono::high_resolution_clock::now();

  int n = DAG->node.size(), m = _seq.size();
  assert(_seq[0] == para->m - 1);
  const std::vector<node_t>& node = DAG->node;const std::vector<int>& rank = DAG->rank;
  std::vector<int> mat = para->mat;
  int match = para->match, mismatch = para->mismatch, para_m = para->m, e1 = para->gap_ext1, o1 = para->gap_open1 + e1;
  reg MATCH = _mm256_set1_epi32(match), MISMATCH = _mm256_set1_epi32(mismatch), E1 = _mm256_set1_epi32(e1), O1 = _mm256_set1_epi32(o1), Neg_inf = _mm256_set1_epi32(NEG_INF);
  // std::cout << mat << " " << mis << " " << o1 << " " << e1 << "\n";
  int reg_size = SIMDTotalBytes / sizeof(int); // number of int in reg
  int col_size = (m + 1 + reg_size - 1) / reg_size * reg_size;
  const size_t mtx_size = n * col_size;
  assert(col_size % reg_size == 0);
  int block_num = col_size / reg_size; // block_num
  std::string seq;seq.resize(col_size);
  for (int i = 0; i < col_size; i++) {
    // if (i < m) seq[(i % block_size) * reg_size + (i / block_size)] = _seq[i];
    // else seq[(i % block_size) * reg_size + (i / block_size)] = char26_table['N'];
    seq[(i % block_num) * reg_size + (i / block_num)] = i < m ? _seq[i] : char26_table['N'];
  }
  // std::cerr << seq.size() << " " << col_size << "\n";
  void* buff = nullptr;
  // ===== 2. 对齐内存分配 =====
  // if ((ret = posix_memalign(&buff, SIMDTotalBytes, 3 * mtx_size * sizeof(int))) != 0) {
  // std::cerr << 3 * mtx_size << "\n";
  // return std::vector<res_t>();
  // std::cerr << buff << "\n";
  if (mpool != nullptr) mpool->alloc_aligned(&buff, SIMDTotalBytes, 3 * mtx_size * sizeof(int));
  else alloc_aligned(&buff, SIMDTotalBytes, 3 * mtx_size * sizeof(int)); // malloc
  int* M = (int*)buff;
  int* D = M + mtx_size;
  int* I = D + mtx_size;
  memset(M, 0xc0, col_size * sizeof(int));//M[0] = NEG_INF
  memset(D, 0xc0, col_size * sizeof(int));//D[0] = NEG_INF
  // memset(buff, 0xc0, 3 * mtx_size * sizeof(int)); // 初始化为-INF 0xc0c0c0c0
  M[0 * col_size + 0] = 0;
  auto end1 = std::chrono::high_resolution_clock::now();
  total_part1 += std::chrono::duration_cast<std::chrono::microseconds>(end1 - start1);
  int tot = 0;
  int offset = (block_num - 1) * reg_size;
  for (int i = 1; i < n; i++) {
    const node_t& cur = node[rank[i]];
    int* M_i = M + i * col_size;
    int* D_i = D + i * col_size;
    tot += cur.in.size();
    auto start2 = std::chrono::high_resolution_clock::now();
    for (int k = 0; k < cur.in.size(); k++) {
      int p = node[cur.in[k]].rank; // rank
      const node_t& pre = node[cur.in[k]];
      reg PRE_BASE = _mm256_set1_epi32(pre.base);
      // j == 0;
      reg Mij = k == 0 ? Neg_inf : _mm256_load_si256((reg*)(M_i));
      reg SEQ; // seq[j - 1]
      // 特殊处理 j=0 的情况
      alignas(32) int tmp[8] = { char26_table['N'], seq[offset + 0], seq[offset + 1], seq[offset + 2],seq[offset + 3], seq[offset + 4], seq[offset + 5], seq[offset + 6] };
      SEQ = _mm256_load_si256((reg*)tmp);

      int* M_p = M + p * col_size;
      // M[i][j] = std::max(M[i][j], M[p][j - 1] + mat[pre.base * para_m + seq[j - 1]]); // qsource
      reg TMP = _mm256_setr_epi32(NEG_INF, M_p[offset + 0], M_p[offset + 1], M_p[offset + 2], M_p[offset + 3], M_p[offset + 4], M_p[offset + 5], M_p[offset + 6]); //M[p][j - 1]
      reg mask = _mm256_cmpeq_epi32(PRE_BASE, SEQ);
      // pre.base == seq[j - 1] ? match : mismatch
      TMP = _mm256_add_epi32(TMP, _mm256_blendv_epi8(MISMATCH, MATCH, mask));
      Mij = _mm256_max_epi32(Mij, TMP); //std::max(M[i][j], M[p][j - 1] + mat[pre.base * para_m + seq[j - 1]]); 
      _mm256_store_si256((reg*)(M_i), Mij); // M[i][j] = max

      for (int j = reg_size; j < col_size; j += reg_size) { // SIMD Msource
        Mij = k == 0 ? Neg_inf : _mm256_load_si256((reg*)(M_i + j));
        SEQ; // seq[j - 1]
        // 直接加载连续内存 (高效)
        const __m128i seq_chunk = _mm_loadu_si128((const __m128i*)(seq.c_str() + j - reg_size));
        SEQ = _mm256_cvtepi8_epi32(seq_chunk);  // 符号扩展到32位

        int* M_p = M + p * col_size;
        // M[i][j] = std::max(M[i][j], M[p][j - 1] + mat[pre.base * para_m + seq[j - 1]]); // qsource
        TMP = _mm256_load_si256((reg*)(M_p + j - reg_size)); //M[p][j - 1]

        mask = _mm256_cmpeq_epi32(PRE_BASE, SEQ);
        // pre.base == seq[j - 1] ? match : mismatch
        TMP = _mm256_add_epi32(TMP, _mm256_blendv_epi8(MISMATCH, MATCH, mask));
        Mij = _mm256_max_epi32(Mij, TMP); //std::max(M[i][j], M[p][j - 1] + mat[pre.base * para_m + seq[j - 1]]); 
        _mm256_store_si256((reg*)(M_i + j), Mij); // M[i][j] = max
      }
      for (int j = 0; j < col_size; j += reg_size) { // SIMD Dsource
        reg Dij = k == 0 ? Neg_inf : _mm256_load_si256((reg*)(D_i + j));
        int* M_p = M + p * col_size;
        int* D_p = D + p * col_size;
        // D[i][j] = std::max({ D[i][j], D[p][j] + e1, M[p][j] + o1 });    // dsource
        TMP = _mm256_load_si256((reg*)(D_p + j));
        Dij = _mm256_max_epi32(Dij, _mm256_add_epi32(TMP, E1)); // max(D[i][j], D[p][j] + e1)
        TMP = _mm256_load_si256((reg*)(M_p + j));
        Dij = _mm256_max_epi32(Dij, _mm256_add_epi32(TMP, O1)); // max(D[i][j], M[p][j] + o1)
        _mm256_store_si256((reg*)(D_i + j), Dij); // D[i][j] = max
      }
    }

    auto end2 = std::chrono::high_resolution_clock::now();
    total_part2 += std::chrono::duration_cast<std::chrono::microseconds>(end2 - start2);
    auto start3 = std::chrono::high_resolution_clock::now();
    // SIMD fsource
    // I[i][0] = NEG_INF;
    int* I_i = I + i * col_size;
    I_i[0] = NEG_INF;
    for (int j = 1; j < reg_size; j++) {
      // I[i][j] = std::max(I[i][j - 1] + e1, M[i][j - 1] + o1)
      I_i[j] = std::max(I_i[j - 1] + block_num * e1, M_i[j - 1 + offset] + o1); // isource
      M_i[j] = std::max({ M_i[j], D_i[j], I_i[j] });  // three source
    }
    for (int j = reg_size; j < col_size; j += reg_size) {
      // I[i][j] = std::max(I[i][j - 1] + e1, M[i][j - 1] + o1)
      reg TMP1 = _mm256_load_si256((reg*)(I_i + j - reg_size));
      TMP1 = _mm256_add_epi32(TMP1, E1);
      reg TMP2 = _mm256_load_si256((reg*)(M_i + j - reg_size));
      TMP2 = _mm256_add_epi32(TMP2, O1);
      _mm256_store_si256((reg*)(I_i + j), _mm256_max_epi32(TMP1, TMP2)); // I[i][j] = max

      // M[i][j] = std::max({ M[i][j], D[i][j], I[i][j] });  // three source
      reg Mij = _mm256_load_si256((reg*)(M_i + j));
      TMP1 = _mm256_load_si256((reg*)(D_i + j));
      TMP2 = _mm256_load_si256((reg*)(I_i + j));
      Mij = _mm256_max_epi32(Mij, _mm256_max_epi32(TMP1, TMP2));
      _mm256_store_si256((reg*)(M_i + j), Mij); // M[i][j] = max
    }
    // active f loop
    I_i[1] = std::max(I_i[offset] + e1, M_i[offset] + o1); // isource
    M_i[1] = std::max({ M_i[1], D_i[1], I_i[1] });  // three source
    for (int j = 2; j < reg_size; j++) {
      // I[i][j] = std::max(I[i][j - 1] + e1, M[i][j - 1] + o1)
      I_i[j] = std::max(I_i[j], I_i[j - 1] + block_num * e1); // isource
      M_i[j] = std::max(M_i[j], I_i[j]);  // three source
    }
    for (int j = reg_size; j < col_size; j += reg_size) {
      // I[i][j] = std::max(I[i][j], I[i][j - 1] + e1)
      reg Iij = _mm256_load_si256((reg*)(I_i + j));
      reg TMP1 = _mm256_load_si256((reg*)(I_i + j - reg_size)); // I[i][j - 1]
      TMP1 = _mm256_add_epi32(TMP1, E1); //I[i][j - 1] + e1
      _mm256_store_si256((reg*)(I_i + j), _mm256_max_epi32(Iij, TMP1)); // I[i][j] = max

      // M[i][j] = std::max(M[i][j], I[i][j]);  // three source
      reg Mij = _mm256_load_si256((reg*)(M_i + j));
      TMP1 = _mm256_load_si256((reg*)(I_i + j));
      _mm256_store_si256((reg*)(M_i + j), _mm256_max_epi32(Mij, TMP1)); // M[i][j] = max
    }
    auto end3 = std::chrono::high_resolution_clock::now();
    total_part3 += std::chrono::duration_cast<std::chrono::microseconds>(end3 - start3);
  }
  // std::cerr << "avg pre size:" << tot / n << "\n";
  // std::cerr << "部分1总耗时: " << total_part1.count() / 1000 << " ms\n";
  // std::cerr << "部分2总耗时: " << total_part2.count() / 1000 << " ms\n";
  // std::cerr << "部分3总耗时: " << total_part3.count() / 1000 << " ms\n";
  // std::cout << M[(n - 1) * col_size + m] << "\n";
  // std::cout << "M:\n";
  // for (int i = 0; i < n; i++) {
  //   for (int j = 0; j < col_size; j++) {
  //     std::cout << M[i * col_size + j] << " \n"[j + 1 == col_size];
  //   }
  // }
  // std::cout << "D:\n";
  // for (int i = 0; i < n; i++) {
  //   for (int j = 0; j < col_size; j++) {
  //     std::cout << D[i * col_size + j] << " \n"[j + 1 == col_size];
  //   }
  // }
  // std::cout << "I:\n";
  // for (int i = 0; i < n; i++) {
  //   for (int j = 0; j < col_size; j++) {
  //     std::cout << I[i * col_size + j] << " \n"[j + 1 == col_size];
  //   }
  // }
  // std::cerr << "finish dp" << "\n";
  int i = n - 1, j = (m % block_num) * reg_size + (m / block_num); // j = m
  // std::cerr << M[i * col_size + j] << "\n";
  std::vector<res_t> res;
  // return res;
  char op = 'M';
  // std::cerr << "finsh align" << "\n";
  while (i > 0 || j > 0) {
    // std::cerr << "op: " << op << "\n";
    // std::cerr << "i: " << i << "\n";
    // std::cerr << "j: " << j << "\n";
    int* M_i = M + i * col_size;
    int* D_i = D + i * col_size;
    int* I_i = I + i * col_size;
    int hit = 0;
    // dsource
    // std::cerr << op << " " << i << " " << j << "\n";
    if (op == 'M') {
      if (M_i[j] == D_i[j]) op = 'D';
      if (M_i[j] == I_i[j]) op = 'I';
    }
    if (op == 'M') {
      const node_t& cur = node[rank[i]];
      for (int k = 0; k < cur.in.size(); k++) {
        int p = node[cur.in[k]].rank; // rank
        const node_t& pre = node[cur.in[k]];
        int* M_p = M + p * col_size;
        if (j - reg_size >= 0) {
          if (pre.base == seq[j - reg_size] && M_i[j] == M_p[j - reg_size] + match) {
            // std::cout << "M";
            res.emplace_back(res_t(pre.id, pre.base));
            i = p, j -= reg_size;
            hit = 1;
            break;
          }
        }
        else {
          if (pre.base == seq[j - 1 + offset] && M_i[j] == M_p[j - 1 + offset] + match) {
            // std::cout << "M";
            res.emplace_back(res_t(pre.id, pre.base));
            i = p, j = j - 1 + offset;
            hit = 1;
            break;
          }
        }
      }

      for (int k = 0; !hit && k < cur.in.size(); k++) {
        int p = node[cur.in[k]].rank; // rank
        const node_t& pre = node[cur.in[k]];
        int* M_p = M + p * col_size;
        if (j - reg_size >= 0) {
          if (pre.base != seq[j - reg_size] && M_i[j] == M_p[j - reg_size] + mismatch) {
            const node_t& par = node[pre.par_id]; // dsu.find par
            if (par.aligned_node[seq[j - reg_size]] != -1) {
              res.emplace_back(res_t(par.aligned_node[seq[j - reg_size]], seq[j - reg_size]));
            }
            else res.emplace_back(res_t(-1, seq[j - reg_size], pre.par_id));
            i = p, j -= reg_size;
            hit = 1;
            break;
          }
        }
        else {
          if (pre.base != seq[j - 1 + offset] && M_i[j] == M_p[j - 1 + offset] + mismatch) {
            // std::cout << "X";
            const node_t& par = node[pre.par_id]; // dsu.find par
            if (par.aligned_node[seq[j - 1 + offset]] != -1) {
              res.emplace_back(res_t(par.aligned_node[seq[j - 1 + offset]], seq[j - 1 + offset]));
            }
            else res.emplace_back(res_t(-1, seq[j - 1 + offset], pre.par_id));
            i = p, j = j - 1 + offset;
            hit = 1;
            break;
          }
        }
      }
      if (hit) continue;
    }
    if (op == 'I') {
      if (j - reg_size >= 0) {
        // std::cout << "I";
        res.emplace_back(res_t(-1, seq[j - reg_size]));
        if (I_i[j] == I_i[j - reg_size] + e1) {
          j -= reg_size;
          continue;
        }
        if (I_i[j] == M_i[j - reg_size] + o1) {
          op = 'M';
          j -= reg_size;
          continue;
        }
      }
      else {
        res.emplace_back(res_t(-1, seq[j - 1 + offset]));
        if (I_i[j] == I_i[j - 1 + offset] + e1) {
          j = j - 1 + offset;
          continue;
        }
        if (I_i[j] == M_i[j - 1 + offset] + o1) {
          op = 'M';
          j = j - 1 + offset;
          continue;
        }
      }
    }
    if (op == 'D') {
      const node_t& cur = node[rank[i]];
      for (int k = 0; k < cur.in.size(); k++) {
        int p = node[cur.in[k]].rank; // rank
        const node_t& pre = node[cur.in[k]];
        int* D_p = D + p * col_size;
        int* M_p = M + p * col_size;
        if (D_i[j] == D_p[j] + e1) {
          i = p;
          break;
        }
      }
      for (int k = 0; k < cur.in.size(); k++) {
        int p = node[cur.in[k]].rank; // rank
        const node_t& pre = node[cur.in[k]];
        int* D_p = D + p * col_size;
        int* M_p = M + p * col_size;
        if (D_i[j] == M_p[j] + o1) {
          op = 'M';
          i = p;
          break;
        }
      }
      continue;
    }
  }
  // std::cerr << "finish POA:" << "\n";
  if (mpool == nullptr) free_aligned(buff);
  return res;
}

std::vector<res_t> POA_SIMD(para_t* para, const graph& DAG, const std::string& _seq, aligned_buff_t* mpool) {
  std::chrono::microseconds total_part1(0);
  std::chrono::microseconds total_part2(0);
  std::chrono::microseconds total_part3(0);
  auto start1 = std::chrono::high_resolution_clock::now();

  int n = DAG.node.size(), m = _seq.size();
  std::string seq = _seq;
  assert(seq[0] == para->m - 1);
  const std::vector<node_t>& node = DAG.node;const std::vector<int>& rank = DAG.rank;
  std::vector<int> mat = para->mat;
  int match = para->match, mismatch = para->mismatch, para_m = para->m, e1 = para->gap_ext1, o1 = para->gap_open1 + e1;
  reg MATCH = _mm256_set1_epi32(match), MISMATCH = _mm256_set1_epi32(mismatch), E1 = _mm256_set1_epi32(e1), O1 = _mm256_set1_epi32(o1), Neg_inf = _mm256_set1_epi32(NEG_INF);
  // std::cout << mat << " " << mis << " " << o1 << " " << e1 << "\n";
  int reg_size = SIMDTotalBytes / sizeof(int); // number of int in reg
  int col_size = (m + 1 + reg_size - 1) / reg_size * reg_size;
  seq.append(col_size - m, char26_table['N']);
  // std::cerr << seq.size() << " " << col_size << "\n";
  assert(col_size % reg_size == 0);
  const size_t mtx_size = n * col_size;
  size_t max_in_size = 0;
  for (int i = 0; i < n; i++) {
    const node_t& cur = node[rank[i]];
    max_in_size = std::max(max_in_size, cur.in.size());
  }
  void* buff = nullptr;
  // ===== 2. 对齐内存分配 =====
  // if ((ret = posix_memalign(&buff, SIMDTotalBytes, 3 * mtx_size * sizeof(int))) != 0) {
  std::cerr << 3 * mtx_size << "\n";
  // return std::vector<res_t>();
  // std::cerr << buff << "\n";
  if (mpool != nullptr) mpool->alloc_aligned(&buff, SIMDTotalBytes, 3 * mtx_size * sizeof(int));
  else alloc_aligned(&buff, SIMDTotalBytes, 3 * mtx_size * sizeof(int)); // malloc
  int* M = (int*)buff;
  int* D = M + mtx_size;
  int* I = D + mtx_size;
  memset(M, 0xc0, col_size * sizeof(int));//M[0] = NEG_INF
  memset(D, 0xc0, col_size * sizeof(int));//D[0] = NEG_INF
  // memset(buff, 0xc0, 3 * mtx_size * sizeof(int)); // 初始化为-INF 0xc0c0c0c0
  M[0 * col_size + 0] = 0;
  auto end1 = std::chrono::high_resolution_clock::now();
  total_part1 += std::chrono::duration_cast<std::chrono::microseconds>(end1 - start1);
  int tot = 0;
  std::vector<int> p(max_in_size), prev(max_in_size);  // 自动初始化值为0 
  reg* PRE_BASE = nullptr;
  alloc_aligned((void**)&PRE_BASE, SIMDTotalBytes, max_in_size * sizeof(reg));
  for (int i = 1; i < n; i++) {
    const node_t& cur = node[rank[i]];
    int* M_i = M + i * col_size;
    int* D_i = D + i * col_size;
    tot += cur.in.size();
    auto start2 = std::chrono::high_resolution_clock::now();
    char prech = char26_table['N'];
    for (int k = 0; k < cur.in.size(); k++) {
      p[k] = node[cur.in[k]].rank; // rank
      prev[k] = NEG_INF;
      const node_t& pre = node[cur.in[k]];
      PRE_BASE[k] = _mm256_set1_epi32(pre.base);
    }
    for (int j = 0; j < col_size; j += reg_size) { // SIMD
      reg Mij = Neg_inf;
      reg Dij = Neg_inf;
      reg SEQ;
      if (j == 0) {
        // 特殊处理 j=0 的情况
        alignas(32) int tmp[8] = {
            prech,
            seq[0], seq[1], seq[2],
            seq[3], seq[4], seq[5], seq[6]
        };
        SEQ = _mm256_load_si256((reg*)tmp);
      }
      else {
        // 直接加载连续内存 (高效)
        const __m128i seq_chunk = _mm_loadu_si128((const __m128i*)(seq.c_str() + j - 1));
        SEQ = _mm256_cvtepi8_epi32(seq_chunk);  // 符号扩展到32位
      }
      for (int k = 0; k < cur.in.size(); k++) {
        int* M_p = M + p[k] * col_size;
        int* D_p = D + p[k] * col_size;
        // M[i][j] = std::max(M[i][j], M[p][j - 1] + mat[pre.base * para_m + seq[j - 1]]); // qsource
    //j - 1;
        reg TMP = j == 0 ? _mm256_setr_epi32(prev[k], M_p[j], M_p[j + 1], M_p[j + 2], M_p[j + 3], M_p[j + 4], M_p[j + 5], M_p[j + 6]) : _mm256_loadu_si256((reg*)(M_p + j - 1)); //M[p][j - 1]

        reg mask = _mm256_cmpeq_epi32(PRE_BASE[k], SEQ);
        // pre.base == seq[j - 1] ? match : mismatch
        TMP = _mm256_add_epi32(TMP, _mm256_blendv_epi8(MISMATCH, MATCH, mask));
        Mij = _mm256_max_epi32(Mij, TMP); //std::max(M[i][j], M[p][j - 1] + mat[pre.base * para_m + seq[j - 1]]); 

        // D[i][j] = std::max({ D[i][j], D[p][j] + e1, M[p][j] + o1 });    // dsource
        TMP = _mm256_load_si256((reg*)(D_p + j));
        Dij = _mm256_max_epi32(Dij, _mm256_add_epi32(TMP, E1)); // max(D[i][j], D[p][j] + e1)
        TMP = _mm256_load_si256((reg*)(M_p + j));
        Dij = _mm256_max_epi32(Dij, _mm256_add_epi32(TMP, O1)); // max(D[i][j], M[p][j] + o1)

        prev[k] = M_p[j + reg_size - 1]; // pre = M[p][j +reg_size - 1]
      }
      prech = seq[j + reg_size - 1];// prech =seq[j + reg_size - 1]
      _mm256_store_si256((reg*)(M_i + j), Mij); // M[i][j] = max
      _mm256_store_si256((reg*)(D_i + j), Dij); // D[i][j] = max
    }
    // if (i % 100 == 0) std::cout << i << "\n";

    auto end2 = std::chrono::high_resolution_clock::now();
    total_part2 += std::chrono::duration_cast<std::chrono::microseconds>(end2 - start2);
    auto start3 = std::chrono::high_resolution_clock::now();
    // I[i][0] = NEG_INF;
    int* I_i = I + i * col_size;
    I_i[0] = NEG_INF;
    for (int j = 1; j <= m; j++) {
      I_i[j] = std::max(I_i[j - 1] + e1, M_i[j - 1] + o1); // isource
      M_i[j] = std::max({ M_i[j], D_i[j], I_i[j] });  // three source
    }
    auto end3 = std::chrono::high_resolution_clock::now();
    total_part3 += std::chrono::duration_cast<std::chrono::microseconds>(end3 - start3);
  }
  std::cerr << "avg pre size:" << tot / n << "\n";
  std::cerr << "部分1总耗时: " << total_part1.count() / 1000 << " ms\n";
  std::cerr << "部分2总耗时: " << total_part2.count() / 1000 << " ms\n";
  std::cerr << "部分3总耗时: " << total_part3.count() / 1000 << " ms\n";
  std::cerr << M[(n - 1) * col_size + m] << "\n";
  std::cout << "M:\n";
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < col_size; j++) {
      std::cout << M[i * col_size + j] << " \n"[j + 1 == col_size];
    }
  }
  std::cout << "D:\n";
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < col_size; j++) {
      std::cout << D[i * col_size + j] << " \n"[j + 1 == col_size];
    }
  }
  std::cout << "I:\n";
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < col_size; j++) {
      std::cout << I[i * col_size + j] << " \n"[j + 1 == col_size];
    }
  }
  // std::cerr << "finish" << "\n";
  int i = n - 1, j = m;
  std::vector<res_t> res;
  char op = 'M';
  // std::cerr << "finsh align" << "\n";
  while (i > 0 || j > 0) {
    // dsource
    std::cerr << op << " " << i << " " << j << "\n";
    if (op == 'M') {
      if (M[i * col_size + j] == D[i * col_size + j]) op = 'D';
      if (M[i * col_size + j] == I[i * col_size + j]) op = 'I';
    }
    // std::cerr << op << " " << i << " " << j << "\n";
    if (op == 'M') {
      const node_t& cur = node[rank[i]];
      for (int k = 0; k < cur.in.size(); k++) {
        int p = node[cur.in[k]].rank; // rank
        const node_t& pre = node[cur.in[k]];
        if (j - 1 >= 0 && M[i * col_size + j] == M[p * col_size + j - 1] + (pre.base == seq[j - 1] ? match : mismatch)) {
          if (pre.base == seq[j - 1]) {
            // std::cout << "M";
            res.emplace_back(res_t(pre.id, pre.base));
          }
          else {
            // std::cout << "X";
            const node_t& par = node[pre.par_id]; // dsu.find par
            if (par.aligned_node[seq[j - 1]] != -1) {
              res.emplace_back(res_t(par.aligned_node[seq[j - 1]], seq[j - 1]));
            }
            else res.emplace_back(res_t(-1, seq[j - 1], pre.par_id));
          }
          i = p, j--;
          break;
        }
      }
      continue;
    }
    if (op == 'I') {
      if (j - 1 >= 0) {
        // std::cout << "I";
        res.emplace_back(res_t(-1, seq[j - 1]));
        if (I[i * col_size + j] == I[i * col_size + j - 1] + e1) {
          j--;
          continue;
        }
        if (I[i * col_size + j] == M[i * col_size + j - 1] + o1) {
          op = 'M';
          j--;
          continue;
        }
      }
    }
    if (op == 'D') {
      const node_t& cur = node[rank[i]];
      // std::cerr << M[i][j] << " " << D[i][j] << "\n";
      for (int k = 0; k < cur.in.size(); k++) {
        int p = node[cur.in[k]].rank; // rank
        const node_t& pre = node[cur.in[k]];
        if (D[i * col_size + j] == D[p * col_size + j] + e1) {
          i = p;
          break;
        }
        if (D[i * col_size + j] == M[p * col_size + j] + o1) {
          op = 'M';
          i = p;
          break;
        }
      }
      continue;
    }
  }
  // std::cout << "sorce:" << M[n - 1][m] << "\n";
  free_aligned(PRE_BASE);
  if (mpool == nullptr) free_aligned(buff);
  return res;
}