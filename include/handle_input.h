#ifndef HANDLE_INPUT_H
#define HANDLE_INPUT_H
#include <vector>
#include "sequence.h"
#include "parameter.h"
#include "graph.h"
void readFile(std::vector<seq_t>& seqs, const char* path);
void initPara(para_t* para);
std::vector<seq_t> read_gfa(para_t* para, graph* DAG, const char* path);
int gfa_parse_S(para_t* para, graph* DAG, seg_seq_t* segs, char* s);
int gfa_parse_W(para_t* para, graph* DAG, std::vector<seq_t>& seqs, seg_seq_t* segs, char* s);
int gfa_parse_P(para_t* para, graph* DAG, std::vector<seq_t>& seqs, seg_seq_t* segs, char* s);
#endif