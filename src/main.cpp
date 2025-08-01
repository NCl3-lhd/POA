#include "handle_input.h"
#include "cxxopts.hpp"
#include <iostream>
#include <zlib.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)
int main(int argc, char** argv) {

  // handle arg
  cxxopts::Options options("POA", "A multiple sequence alignment tool");

  options.add_options()
    ("i,input", "input path", cxxopts::value<std::string>())
    ("v,verbose", "detailed output mode", cxxopts::value<bool>()->default_value("false"))
    ("h,help", "Print usage")
    ;
  std::string path;
  try {
    auto result = options.parse(argc, argv);
    if (result.count("help"))
    {
      std::cout << options.help() << std::endl;
      return 0;
    }
    if (result.count("input"))
      path = result["input"].as<std::string>();
    std::cout << path << "\n";
  }
  catch (const cxxopts::exceptions::exception& e)
  {
    std::cout << "error parsing options: " << e.what() << std::endl;
    return 1;
  }

  //handle read file
  gzFile fp; // 文件指针
  kseq_t* seq; // 序列结构体
  int l; // 用于存储读取序列的长度

  // 检查命令行参数，确保提供了输入文件
  if (path.empty()) {
    fprintf(stderr, "Usage: %s <in.fasta>\n", argv[0]);
    return 1;
  }

  // 打开 Gzip 压缩的 FASTA 文件
  fp = gzopen(path.c_str(), "r");
  if (fp == NULL) {
    fprintf(stderr, "Error opening file %s\n", argv[1]);
    return 1;
  }

  // 初始化 kseq 结构体
  seq = kseq_init(fp);

  // 逐行读取序列
  while ((l = kseq_read(seq)) >= 0) {
    printf("name: %s\n", seq->name.s); // 打印序列名称
    if (seq->comment.l) printf("comment: %s\n", seq->comment.s); // 打印序列注释（如果存在）
    printf("seq: %s\n", seq->seq.s); // 打印序列
    if (seq->qual.l) printf("qual: %s\n", seq->qual.s); // 打印质量分数（如果存在）
  }
  // 释放 kseq 结构体
  kseq_destroy(seq);

  // 关闭文件
  gzclose(fp);

  return 0;
}