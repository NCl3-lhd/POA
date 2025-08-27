#ifndef MINIMIZER_H
#define MINIMIZER_H
#include "parameter.h"
#include "sequence.h"
#include <cstdint>

typedef struct { uint64_t x, y; } mm128_t;
typedef struct { size_t n, m; mm128_t* a; } mm128_v;
void mm_sketch(void* km, const char* str, int len, int w, int k, uint32_t rid, int is_hpc, mm128_v* p);
struct minimizer_t {
  void* km;
  size_t seqs_size;
  int maxl;
  std::vector<int> len;
  int* mm_h;
  mm128_v mm_v;
  std::vector<int> max_sim; // is seq_id;
  minimizer_t(para_t* para, const std::vector<seq_t>& seqs);
  ~minimizer_t();
  void init(para_t* para, const std::vector<seq_t>& seqs);
  int collect_mm(void* km, const std::vector<seq_t>& seqs, para_t* para);
  std::vector<int> get_guide_tree(para_t* para);
};
// std::vector<sequence> readFile(const char* path);
#endif