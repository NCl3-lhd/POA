#include "align.h"
#include "cassert"
#include <algorithm>
constexpr int INF = 100; // 0x3f3f3f3f
void POA(para_t* para, const graph& DAG, const std::string& seq) {
  int n = DAG.node.size(), m = seq.size();
  assert(seq[0] == para->m - 1);
  const std::vector<node_t>& node = DAG.node;const std::vector<int>& rank = DAG.rank;
  std::vector<int> mat = para->mat; int para_m = para->m, o1 = para->gap_open1, e1 = para->gap_ext1;
  // std::cout << mat << " " << mis << " " << o1 << " " << e1 << "\n";
  std::vector<std::vector<int>> M(n, std::vector<int>(m + 1, -INF)), D(M), I(M);
  M[0][0] = 0;
  std::cout << n << " " << m << "\n";
  for (int i = 0; i < n; i++) {
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
  std::cout << M[n - 1][m] << "\n";
}