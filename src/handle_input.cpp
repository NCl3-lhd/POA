#include <zlib.h>
#include <stdexcept> 
#include "handle_input.h"
#include "kseq.h"
#include "sequence.h" 

KSEQ_INIT(gzFile, gzread)

// handle input file
std::vector<seq_t> readFile(const char* path) {
  //handle read file
  gzFile fp; // 文件指针
  kseq_t* seq; // 序列结构体
  int l; // 用于存储读取序列的长度
  // 打开 Gzip 压缩的 FASTA 文件
  fp = gzopen(path, "r");
  if (fp == nullptr) {
    throw std::runtime_error("Error opening file " + std::string(path));
  }

  // 初始化 kseq 结构体
  seq = kseq_init(fp);

  // 逐行读取序列
  std::vector<seq_t> seqs;
  while ((l = kseq_read(seq)) >= 0) {
    // printf("name: %s\n", seq->name.s); // 打印序列名称
    seq_t tseq;
    tseq.name = seq->name.s;
    if (seq->comment.l) tseq.comment = seq->comment.s;// 序列注释
    tseq.seq = seq->seq.s;
    if (seq->qual.l)  tseq.qual = seq->qual.s;  // 质量分数
    seqs.emplace_back(tseq);
  }
  // 释放 kseq 结构体
  kseq_destroy(seq);
  // 关闭文件
  gzclose(fp);
  return seqs;
}
// handle input arg
void initPara(para_t* para) {
  //handle incorrect para
  if (para->match < 0) para->m *= -1;
  if (para->mismatch > 0) para->mismatch *= -1;
  if (para->gap_open1 > 0) para->gap_open1 *= -1;
  if (para->gap_ext1 > 0) para->gap_ext1 *= -1;
  para->m = 5; // default m = Nucleotide num
  if (para->mat_fp.empty()) {
    int m = para->m;
    para->mat.resize(m * m, 0); //like HOXD70.mtx
    for (int i = 0; i < m - 1; i++) {
      for (int j = 0; j < m - 1; j++) {
        para->mat[i * m + j] = i == j ? para->match : para->mismatch;
      }
    }
  }
  else { // 待支持 

  }
  if (para->m > 5) { // for aa sequence
    for (int i = 0; i < 256; ++i) {
      char26_table[i] = aa26_table[i];
      char256_table[i] = aa256_table[i];
    }
  }
  else {
    for (int i = 0; i < 256; ++i) {
      char26_table[i] = nt4_table[i];
      char256_table[i] = nt256_table[i];
    }
  }

}