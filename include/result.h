#ifndef RESULT_H
#define RESULT_H
struct res_t {
  // 成员声明
  int from; // if from < 0, new node
  unsigned char base;
  int aligned_id = -1;
  res_t(int _from, unsigned char _base, int _aligned_id = -1) {
    from = _from;
    base = _base;
    aligned_id = _aligned_id;
  }
};
#endif