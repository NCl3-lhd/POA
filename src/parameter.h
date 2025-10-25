#ifndef PARAMETER_H
#define PARAMETER_H
#include <vector>
#include <string>
struct para_t {
  // 成员声明
  int m;std::vector<int> mat;std::string mat_fp;// score matrix  m is alphabet size
  int match, mismatch, gap_open1, gap_open2, gap_ext1, gap_ext2; int inf_min;
  int b, f;// adband
  int k, w;// k-mer
  bool enable_seeding;
  bool progressive_poa;
  int result; // result format 0 consensus,1 rc-msa
  int verbose;  //verbose level (0-2). 0: none, 1: information, 2: debug
};
#endif