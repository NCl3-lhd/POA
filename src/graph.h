#ifndef GRAPH_H
#define GRAPH_H
#include <vector>
#include <string>
#include "node.h"
struct graph {
  // 成员声明
  std::vector<node_t> node; //index is id
  std::vector<int> rank; // rank_to_node_id
  void init(const std::string& str, int para_m, int seq_id);
  void add_adj(int from, int to, int seq_id);
  void topsort();
};
// std::vector<sequence> readFile(const char* path);
#endif