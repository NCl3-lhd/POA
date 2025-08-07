#include "align.h"
#include "cassert"
#include <algorithm>
constexpr int INF = 0x3f3f3f3f; // 0x3f3f3f3f
std::vector<res_t> POA(para_t* para, const graph& DAG, const std::string& seq) {
  int n = DAG.node.size(), m = seq.size();
  assert(seq[0] == para->m - 1);
  const std::vector<node_t>& node = DAG.node;const std::vector<int>& rank = DAG.rank;
  std::vector<int> mat = para->mat; int para_m = para->m, e1 = para->gap_ext1, o1 = para->gap_open1 + e1;
  // std::cout << mat << " " << mis << " " << o1 << " " << e1 << "\n";
  std::vector<std::vector<int>> M(n, std::vector<int>(m + 1, -INF)), D(M), I(M);
  M[0][0] = 0;
  // std::cout << "n:" << n << " " << "m:" << m << "\n";
  for (int i = 1; i < n; i++) {
    const node_t& cur = node[rank[i]];
    // if (i % 100 == 0) std::cout << i << "\n";
    for (int j = 0; j <= m; j++) {
      for (int k = 0; k < cur.in.size(); k++) {
        int p = node[cur.in[k]].rank; // rank
        const node_t& pre = node[cur.in[k]];
        if (j - 1 >= 0) M[i][j] = std::max(M[i][j], M[p][j - 1] + mat[pre.base * para_m + seq[j - 1]]); // qsource
        D[i][j] = std::max({ D[i][j], D[p][j] + e1, M[p][j] + o1 });    // dsource
      }
      if (j - 1 >= 0) I[i][j] = std::max(I[i][j - 1] + e1, M[i][j - 1] + o1); // isource
      M[i][j] = std::max({ M[i][j], D[i][j], I[i][j] });  // three source
    }
  }

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
        if (j - 1 >= 0 && M[i][j] == M[p][j - 1] + mat[pre.base * para_m + seq[j - 1]]) {
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
}