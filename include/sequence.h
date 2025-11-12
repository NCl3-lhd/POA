#ifndef SEQUENCE_H
#define SEQUENCE_H
#include <string>
#include "kstring.h"
#include "khash.h"
#include "kvec.h"

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

KHASH_MAP_INIT_STR(str, uint32_t)
typedef struct {
  int n, m;
  kstring_t* seq, * name;
  khash_t(str)* h; // name -> seq_pos in seq
  kvec_t(uint32_t) start_id, end_id;  // id is the node id
} seg_seq_t;

seg_seq_t* seg_seq_init(void);
seg_seq_t* seg_seq_realloc(seg_seq_t* r);
void seg_seq_free(seg_seq_t* s);
#endif