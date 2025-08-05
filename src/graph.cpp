#include "graph.h"
#include "sequence.h"
#include <queue>
extern char char26_table[256];
void graph::add_adj(int seq_id, int from, int to) {
  // std::cout << from << " " << to << "\n";
  node[from].add_out_adj(seq_id, to);
  node[to].add_in_adj(seq_id, from);
}
void graph::topsort() {
  std::queue<int> q;
  std::vector<int> deg(node.size());
  for (int i = 0; i < node.size(); i++) {
    deg[i] = node[i].in.size();
    if (deg[i] == 0) q.push(i);
  }
  rank.clear();
  while (q.size()) {
    int u = q.front();
    q.pop();
    rank.emplace_back(u);
    const node_t& cur = node[u];
    for (int i = 0; i < cur.out.size(); i++) {
      int v = cur.out[i];
      if (--deg[v] == 0) {
        q.push(v);
      }
    }
  }
  for (int i = 0; i < rank.size(); i++) {
    node[rank[i]].rank = i;
  }
}
void graph::init(const std::string& str, int para_m, int seq_id) {
  // 清空现有数据
  node.clear();
  node.emplace_back(node_t(node.size(), char26_table['N'], para_m)); // src 0
  node.emplace_back(node_t(node.size(), char26_table['N'], para_m)); // sink 1
  int pre_id = 0;
  for (int i = 0; i < str.size(); i++) {
    int cur_id = node.size();
    node.emplace_back(node_t(cur_id, char26_table[str[i]], para_m)); // str[i]
    add_adj(seq_id, pre_id, cur_id);
    pre_id = cur_id;
  }
  add_adj(seq_id, pre_id, 1); // 
  topsort();
}
