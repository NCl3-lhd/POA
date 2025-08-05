#include <iostream>
#include "cxxopts.hpp"
#include "handle_input.h"
#include "parameter.h"
#include "graph.h"
#include "align.h"
// #include "alignment"

// extern unsigned char nt4_table[256];
int main(int argc, char** argv) {

  // handle arg
  cxxopts::Options options("POA", "A multiple sequence alignment tool");

  options.add_options()
    ("i,input", "input path", cxxopts::value<std::string>())
    ("m,mat_fp", "match file path", cxxopts::value<std::string>())
    ("M,match", "match sorce", cxxopts::value<int>()->default_value("1"))
    ("X,mismatch", "mismatch sorce", cxxopts::value<int>()->default_value("-2"))
    ("O,gap_open", "gap_open sorce", cxxopts::value<int>()->default_value("-3"))
    ("E,gap_ext", "gap_ext sorce", cxxopts::value<int>()->default_value("-1"))
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
  std::cout << seqs.size() << "\n";
  // handle alignment 
  graph DAG;
  // seqs[0].seq = "GAAAAGATTTACATAATACTTAGAAAATTTTCCTGGTGGGTTAACTCTGAGCTATGATTTTTTAA";
  // seqs[1].seq = "GAATGATTATGATATACTTAGAAAATTTTAAGTTAATGCACTGTTTCCGACAAATGTGATGATTTTTTAATGATT";
  DAG.init(seqs[0].seq, para->m, 0);
  std::cout << seqs[0].seq << "\n" << seqs[1].seq << "\n";
  for (int i = 1; i < 2; i++) {  //seqs.size()
    std::string tseq;
    tseq += char26_table['N'];
    for (int j = 0; j < seqs[i].seq.size(); j++) {
      // std::cout << (int)char26_table[seqs[i].seq[j]] << " ";
      tseq += char26_table[seqs[i].seq[j]];
    }
    POA(para, DAG, tseq);
    // DAG.add_path(ret);
    // DAG.topsort();
  }

  // handle output 


  // delete
  delete para;
  para = nullptr;  // 防止后续误用
  return 0;
}