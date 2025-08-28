#ifndef GRAPH_H
#define GRAPH_H
#include <vector>
#include <string>
#include "node.h"
#include "result.h"
#include "sequence.h"
#include "parameter.h"
struct graph {
  // 成员声明
  std::vector<node_t> node; //index is id
  // std::vector<int> node_h;  //index is ord[id] 
  std::vector<int> rank; // rank_to_node_id
  std::vector<int> hmin, hmax, tmin, tmax;  // index is rank  haid tail
  std::vector<int> hlen, tlen;  // index is rank  haid tail
  void init(para_t* para, int seq_id, const std::string& str);
  void add_adj(int from, int to, int seq_id);
  void add_path(int para_m, int seq_id, const std::vector<res_t>& res, int graph_node_num = 0);
  void topsort(int op, int para_f);
  void output_rc_msa(const std::vector<seq_t>& seqs);
  std::vector<int> calculateR() const;
};
// std::vector<sequence> readFile(const char* path);
#endif