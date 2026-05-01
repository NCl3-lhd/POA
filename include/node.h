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
  // std::vector<int> ids; // restore the ord[rid]  msa
  // std::vector<int> idp; // restore the ord[rid] pos  
  node_t();
  node_t(int _id, unsigned char _base, int m) {
    base = _base;
    id = par_id = _id;
    ind = 0;  
    aligned_node.resize(m, -1);
    aligned_node[base] = id;
  }
  void add_in_adj(int from) { // seq_id == ord
    int ok = 0;
    ind++;
    for (size_t i = 0; i < in.size(); i++) {
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
  }
  void add_out_adj(int to) {
    int ok = 0;
    for (size_t i = 0; i < out.size(); i++) {
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
  }
  // int getPos(int ord) const {
  //   size_t idx = std::lower_bound(ids.begin(), ids.end(), ord) - ids.begin();
  //   return idp[idx];
  //   if (idx < ids.size() && ids[idx] == ord) return idp[idx];
  //   return -1; // INF
  // }
};
// std::vector<sequence> readFile(const char* path);
#endif