#ifndef SEQUENCE_H
#define SEQUENCE_H
#include <string>


extern unsigned char nt4_table[256];
extern char nt256_table[256];
extern unsigned char aa26_table[256];
extern char aa256_table[256];
extern char char26_table[256];
extern char char256_table[256];
struct seq_t {
  // 成员声明
  std::string name, comment, seq, qual;
};
#endif