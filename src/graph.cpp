#include "graph.h"
#include "sequence.h"
#include <queue>
extern char char26_table[256];
void graph::add_adj(int seq_id, int from, int to) {
  // from and to node actually exist
  // std::cout << from << " " << to << "\n";
  if (from == to) return; // self loop is error
  node[from].add_out_adj(seq_id, to);
  node[to].add_in_adj(seq_id, from);
}
void graph::topsort(int op) { // if op == 1, is not normal topsort, the node_id in queue must be the parent
  std::queue<int> q;
  std::vector<int> deg(node.size());
  if (op == 1) {
    for (int i = 0; i < node.size(); i++) {
      if (node[i].par_id != node[i].id) continue;// is not parent
      int sum_deg = 0;
      for (int j = 0; j < node[i].aligned_node.size(); j++) {
        if (node[i].aligned_node[j] < 0) continue;
        sum_deg += node[node[i].aligned_node[j]].in.size();
      }
      deg[i] = sum_deg;
      // std::cout << i << " " << node[i].par_id << " " << sum_deg << "\n";
      if (deg[i] == 0) q.push(i);
    }
    rank.clear();
    while (q.size()) {
      int u = q.front();
      // std::cout << u << "\n";
      q.pop();
      node_t& cur = node[u];
      cur.rank = rank.size();
      rank.emplace_back(u);
      for (int i = 0; i < cur.aligned_node.size(); i++) {
        if (cur.aligned_node[i] < 0) continue;
        const node_t& son = node[cur.aligned_node[i]];
        for (int j = 0; j < son.out.size(); j++) {
          int v = node[son.out[j]].par_id;
          if (--deg[v] == 0) {
            q.push(v);
          }
        }
      }
    }
    return;
  }
  for (int i = 0; i < node.size(); i++) {
    deg[i] = node[i].in.size();
    if (deg[i] == 0) q.push(i);
  }
  rank.clear();
  while (q.size()) {
    int u = q.front();
    q.pop();
    node_t& cur = node[u];
    cur.rank = rank.size();
    rank.emplace_back(u);
    for (int i = 0; i < cur.out.size(); i++) {
      int v = cur.out[i];
      if (--deg[v] == 0) {
        q.push(v);
      }
    }
  }
  hlen.resize(node.size());
  for (int i = 0; i < rank.size(); i++) { // ni topsort id dp
    int u = rank[i];
    const node_t& cur = node[u];
    int wmax = -1, max_pre = -1;
    for (int k = 0; k < cur.in.size(); k++) {
      int v = cur.in[k];
      int pre = node[v].rank;
      if (cur.in_weight[k] > wmax) {
        wmax = cur.in_weight[k];
        max_pre = pre;
      }
    }
    if (max_pre != -1) hlen[cur.rank] = hlen[max_pre] + 1;
    // std::cerr << u << " " << lp[i] << " " << rp[i] << "\n";
  }

  tlen.resize(node.size());
  tlen[node[1].rank] = 1;
  for (int i = int(rank.size()) - 1; i >= 0; i--) { // ni topsort id dp
    int u = rank[i];
    const node_t& cur = node[u];
    int wmax = -1, max_suc = -1;
    for (int k = 0; k < cur.out.size(); k++) {
      int v = cur.out[k];
      int suc = node[v].rank;
      if (cur.out_weight[k] > wmax) {
        wmax = cur.out_weight[k];
        max_suc = suc;
      }
    }
    if (max_suc != -1) tlen[cur.rank] = tlen[max_suc] + 1;

    // std::cerr << u << " " << lp[i] << " " << rp[i] << "\n";
  }
  // lp[node[0].rank] = 0; // src lp = 0
}
void graph::init(int para_m, int seq_id, const std::string& str) {
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
  topsort(0);
}
void graph::add_path(int para_m, int seq_id, const std::vector<res_t>& res, int graph_node_num) {
  int anchored_id = 1; // sink
  // std::cerr << seq_id << " " << res.size() << "\n";
  for (int i = 0; i < res.size(); i++) {
    int cur_id = res[i].from;
    if (cur_id < 0) {
      if (res[i].aligned_id >= 0) { // dsu.merge
        if (node[res[i].aligned_id].aligned_node[res[i].base] == -1) {
          // std::cerr << 1 << "\n";
          cur_id = node.size();
          node.emplace_back(node_t(cur_id, res[i].base, para_m));
          node[res[i].aligned_id].aligned_node[res[i].base] = cur_id;
          node[cur_id].par_id = res[i].aligned_id;
          add_adj(seq_id, cur_id, anchored_id);
        }
        else {
          // std::cerr << 2 << "\n";
          cur_id = node[res[i].aligned_id].aligned_node[res[i].base];
          add_adj(seq_id, cur_id, anchored_id);
        }
      }
      else {
        // std::cerr << 3 << "\n";
        if (graph_node_num > 0) {
          for (int k = 0; k < node[anchored_id].in.size(); k++) {
            int v = node[anchored_id].in[k];
            if (node[v].base == res[i].base && v >= graph_node_num && node[v].par_id == v) {
              cur_id = v;
              break;
            }
          }
          // for (int k = 0; k < node[anchored_id].in.size(); k++) {
          //   int v = node[anchored_id].in[k];
          //   if ((node[v].base ^ 2) == res[i].base && v >= graph_node_num) {
          //     cur_id = v;
          //   }
          // }
        }
        if (cur_id < 0) {
          cur_id = node.size();
          node.emplace_back(node_t(cur_id, res[i].base, para_m));
        }
        add_adj(seq_id, cur_id, anchored_id);
      }
      // std::cout << cur_id << " " << char256_table[res[i].base] << "  " << anchored_id << "\n";
    }
    else {
      // std::cerr << 4 << "\n";
      // std::cerr << node.size() << " " << cur_id << " " << anchored_id << "\n";
      // std::cout << cur_id << " " << char256_table[res[i].base] << "  " << anchored_id << "\n";
      add_adj(seq_id, cur_id, anchored_id);
    }
    anchored_id = cur_id;
  }
  // std::cerr << "finish add path" << "\n";
}

void graph::output_rc_msa(const std::vector<seq_t>& seqs) {
  // std::cerr << seqs.size() << " " << rank.size() << "\n";
  std::vector<std::string> res(seqs.size(), std::string(rank.size() - 2, '-'));
  for (int i = 0; i < node.size(); i++) {
    int rank = node[node[i].par_id].rank;
    for (int id : node[i].ids) {
      if (rank) res[id][rank - 1] = char256_table[node[i].base];
    }
  }
  for (int i = 0; i < seqs.size(); i++) {
    std::cout << ">" << seqs[i].name << " " << seqs[i].comment << "\n";
    std::cout << res[i] << "\n";
  }
}
std::vector<int> graph::calculateR() const {
  std::queue<int> q;
  std::vector<int> deg(node.size());
  std::vector<int> R(node.size());  // index is rank
  for (int i = 0; i < node.size(); i++) {
    R[i] = 0;
    deg[i] = node[i].out.size();
  }
  q.push(1);  // sink
  while (q.size()) {
    int u = q.front();
    q.pop();
    const node_t& cur = node[u];
    if (u == 1) {
      R[cur.rank] = -1;
    }
    else {
      int wmax = -1, max_suc = -1;
      for (int k = 0; k < cur.out.size(); k++) {
        int v = cur.out[k];
        int suc = node[v].rank;
        if (cur.out_weight[k] > wmax) {
          wmax = cur.out_weight[k];
          max_suc = suc;
        }
      }
      R[cur.rank] = R[max_suc] + 1;
    }
    if (u == 0) {
      break;
    }
    for (int k = 0; k < cur.in.size(); k++) {
      int v = node[u].in[k];
      if (--deg[v] == 0) {
        q.push(v);
      }
    }
  }
  return R;
}