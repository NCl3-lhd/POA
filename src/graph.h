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
  void init(para_t* para);
  void init(para_t* para, int seq_id, const std::string& str);
  int add_node(para_t* para, char base);
  void add_adj(int seq_id, int from, int to, int curPos);
  void add_path(int para_m, int seq_id, const std::vector<res_t>& res, int graph_node_num = 0);
  void topsort(const para_t* para, int op);
  void output_rc_msa(para_t* para, const std::vector<int>& rid_to_ord, const std::vector<seq_t>& seqs);
  void output_consensus();
  void build_consensus();
  void output_gfa(const std::vector<int>& rid_to_ord, const std::vector<seq_t>& seqs);
  std::vector<int> calculateR() const;
  bool is_topsorted; //
  std::string cons;
  std::vector<int> cons_pos_to_id;
};
// std::vector<sequence> readFile(const char* path);
#endif