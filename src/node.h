#ifndef NODE_H
#define NODE_H
#include <vector>
struct node {
  // 成员声明
  int rank;
  char base;
  std::vector<int> in, out;
  int isPar;
  std::vector<int> aligned_node;


};
// std::vector<sequence> readFile(const char* path);
#endif