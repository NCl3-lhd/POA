#include <iostream>
#include "cxxopts.hpp"
#include "handle_input.h"
#include "parameter.h"
#include "graph.h"
#include "align.h"
#include <immintrin.h>
#include "mem_alloc_utils.h"
// #include "kband.h"

// extern unsigned char nt4_table[256];
int main(int argc, char** argv) {

  // handle arg
  cxxopts::Options options("POA", "A multiple sequence alignment tool");

  options.add_options()
    ("i,input", "input path", cxxopts::value<std::string>())
    ("m,mat_fp", "match file path", cxxopts::value<std::string>())
    ("M,match", "match sorce", cxxopts::value<int>()->default_value("2"))
    ("X,mismatch", "mismatch sorce", cxxopts::value<int>()->default_value("-4"))
    ("O,gap_open", "gap_open sorce", cxxopts::value<int>()->default_value("-4"))
    ("E,gap_ext", "gap_ext sorce", cxxopts::value<int>()->default_value("-2"))
    ("h,help", "Print usage")
    ;
  std::string path;
  para_t* para = new para_t();
  try {
    auto result = options.parse(argc, argv);
    if (result.count("help"))
    {
      std::cout << options.help() << std::endl;
      return 0;
    }
    if (result.count("input"))
      path = result["input"].as<std::string>();
    if (result.count("mat_fp"))
      para->mat_fp = result["mat_fp"].as<std::string>();
    para->match = result["match"].as<int>();
    para->mismatch = result["mismatch"].as<int>();
    para->gap_open1 = result["gap_open"].as<int>();
    para->gap_ext1 = result["gap_ext"].as<int>();
  }
  catch (const cxxopts::exceptions::exception& e)
  {
    std::cerr << "error parsing options: " << e.what() << std::endl;
    return 1;
  }
  initPara(para);

  // handle input
  std::vector<seq_t> seqs;
  try {
    seqs = readFile(path.c_str());
  }
  catch (const std::exception& e) {
    std::cerr << "error read file: " << e.what() << std::endl;
    return 1;
  }
  std::cerr << seqs.size() << "\n";
  // handle alignment 
  graph DAG;
  DAG.init(para->m, 0, seqs[0].seq);
  // seqs[0].seq = "TTGCCCTT";
  // seqs[1].seq = "CCAATTTT";
  // seqs[2].seq = "TGCT";
  aligned_buff_t mpool;
  for (int i = 1; i < seqs.size(); i++) {  //seqs.size()
    std::cerr << i << "\n";
    std::string tseq;
    tseq += char26_table['N'];
    for (int j = 0; j < seqs[i].seq.size(); j++) {
      tseq += char26_table[seqs[i].seq[j]];
    }
    // std::cerr << "poa" << "\n";
    // POA_SIMD(para, DAG, tseq);
    // std::vector<res_t> res = POA(para, DAG, tseq);
    // std::vector<res_t> res = POA_SIMD(para, DAG, tseq);
    std::vector<res_t> res = POA_SIMD(para, DAG, tseq, &mpool);
    // std::cerr << "add_path" << "\n";
    DAG.add_path(para->m, i, res);
    // std::cerr << "topsort" << "\n";
    DAG.topsort(i + 1 == seqs.size());
    // std::cout << i << " " << DAG.rank.size() << "\n";
  }
  // handle output 
  DAG.output_rc_msa(seqs);


  // std::cout << "correct check" << "\n";
  // std::string s1 = "TTGCCCTT";
  // std::string s2 = "CCAATTTT";
  // std::string s3 = "CCTT";
  // std::string s4 = "TTGCCCAATTTT";
  // std::string alignedS, alignedT;
  // std::cout << PSA_Kband(s1, seqs[2].seq, &alignedS, &alignedT) << "\n";
  // std::cout << alignedS << " " << alignedT << "\n";
  // std::cout << PSA_Kband(s2, seqs[2].seq, nullptr, nullptr) << "\n";
  // std::cout << PSA_Kband(s3, seqs[2].seq, nullptr, nullptr) << "\n";
  // std::cout << PSA_Kband(s4, seqs[2].seq, nullptr, nullptr) << "\n";
  // delete
  delete para;
  para = nullptr;  // 防止后续误用
  return 0;
}