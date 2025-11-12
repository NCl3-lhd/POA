#ifndef NODE_H
#define NODE_H
#include <vector>
#include <iostream>
struct node_t {
  // 成员声明
  int id, rank, par_id;
  unsigned char base;
  int ind; // in_degree
  std::vector<int> in, in_weight, out, out_weight;
  std::vector<int> aligned_node;
  std::vector<int> ids; // restore the ord[rid]
  std::vector<int> idp; // restore the ord[rid] pos 
  std::vector<int> idpos;
  node_t();
  node_t(int _id, unsigned char _base, int m) {
    base = _base;
    id = par_id = _id;
    ind = 0;  
    aligned_node.resize(m, -1);
    aligned_node[base] = id;
  }
  void add_in_adj(int seq_id, int from, int curPos) { // seq_id == ord
    int ok = 0;
    ind++;
    for (int i = 0; i < in.size(); i++) {
      if (from == in[i]) {
        in_weight[i]++;
        ok = 1;
        break;
      }
    }
    if (!ok) {
      in.emplace_back(from);
      in_weight.emplace_back(1);
    }
    if (ids.empty() || ids.back() != seq_id) {
      // std::cerr << seq_id << " " << (int)base << " " << base << " " << curPos << "\n";
      ids.emplace_back(seq_id);
      idp.emplace_back(curPos);
    }
  }
  void add_out_adj(int seq_id, int to) {
    int ok = 0;
    for (int i = 0; i < out.size(); i++) {
      if (to == out[i]) {
        out_weight[i]++;
        ok = 1;
        break;
      }
    }
    if (!ok) {
      out.emplace_back(to);
      out_weight.emplace_back(1);
    }
    // if (seq_id == 2 && to == 14) std::cerr << id << " " << seq_id << "\n";
    if (ids.empty() || ids.back() != seq_id) {
      ids.emplace_back(seq_id);
    }

  }
  int getPos(int ord) const {
    int idx = std::lower_bound(ids.begin(), ids.end(), ord) - ids.begin();
    return idp[idx];
    if (idx < ids.size() && ids[idx] == ord) return idp[idx];
    return -1; // INF
  }
};
// std::vector<sequence> readFile(const char* path);
#endif