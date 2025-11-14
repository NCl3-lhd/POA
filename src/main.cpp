#include <iostream>
#include <algorithm>
#include <numeric>
#include "cxxopts.hpp"
#include "minipoa.h"
// #include "kband.h"

// extern unsigned char nt4_table[256];
int main(int argc, char** argv) {

  // handle arg
  cxxopts::Options options("POA", "A multiple sequence alignment tool");

  options.add_options()
    ("input", "input path", cxxopts::value<std::string>())
    ("i,inc_fp", "incrementally align sequences to an existing graph/MSA", cxxopts::value<std::string>())
    ("m,mat_fp", "match file path", cxxopts::value<std::string>())
    ("M,match", "match sorce", cxxopts::value<int>()->default_value("2"))
    ("X,mismatch", "mismatch sorce", cxxopts::value<int>()->default_value("-4"))
    ("O,gap_open", "gap_open sorce", cxxopts::value<int>()->default_value("-4"))
    ("E,gap_ext", "gap_ext sorce", cxxopts::value<int>()->default_value("-2"))
    ("t,thread", "thread number", cxxopts::value<int>()->default_value("1"))
    ("b,band_b", "band arg", cxxopts::value<int>()->default_value("100"))
    ("f,band_f", "band arg", cxxopts::value<int>()->default_value("40"))
    ("B,ab_band", "adpative band arg", cxxopts::value<bool>()->default_value("false"))
    ("S,seeding", "enable minimizer-based seeding and anchoring", cxxopts::value<bool>()->default_value("false"))
    ("W,poa_w", "the minimum distance between adjacent anchors", cxxopts::value<int>()->default_value("10000"))
    ("k,k_mer", "k_mer lenth", cxxopts::value<int>()->default_value("19")) //19
    ("w,mm_w", "k_mer_window lenth", cxxopts::value<int>()->default_value("10"))
    ("s,sample_num", "sample_num", cxxopts::value<int>()->default_value("50"))
    ("p,progressive_poa", "is progressive_poa", cxxopts::value<bool>()->default_value("false"))
    ("r,result", "result format (0-2). 0: consensus, 1: rc-msa, 2: gfa", cxxopts::value<int>()->default_value("0"))
    ("V,verbose", "verbose level (0-2). 0: none, 1: information, 2: debug [0]\n", cxxopts::value<int>()->default_value("0"))
    ("h,help", "Print usage")
    ;
  options.parse_positional({ "input" });
  int thread, sample_num;
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
    if (result.count("inc_fp"))
      para->inc_fp = result["inc_fp"].as<std::string>();
    para->poa_w = 0;
    if (result.count("poa_w"))
      para->poa_w = result["poa_w"].as<int>();
    para->match = result["match"].as<int>();
    para->mismatch = result["mismatch"].as<int>();
    para->gap_open1 = result["gap_open"].as<int>();
    para->gap_ext1 = result["gap_ext"].as<int>();
    para->b = result["band_b"].as<int>();
    para->f = result["band_f"].as<int>();
    para->ab_band = result["ab_band"].as<bool>();
    para->k = result["k_mer"].as<int>();
    para->mm_w = result["mm_w"].as<int>();
    para->bw = 1000;
    para->progressive_poa = result["progressive_poa"].as<bool>();
    para->enable_seeding = result["seeding"].as<bool>();
    para->thread = result["thread"].as<int>();
    sample_num = result["sample_num"].as<int>();
    if (sample_num < 0) sample_num *= -1;
    para->result = result["result"].as<int>();
    para->verbose = result["verbose"].as<int>();

  }
  catch (const cxxopts::exceptions::exception& e)
  {
    std::cerr << "error parsing options: " << e.what() << std::endl;
    return 1;
  }
  initPara(para);

  // handle input
  std::vector<seq_t> seqs;
  graph* DAG = new graph();
  DAG->init(para);
  if (!para->inc_fp.empty()) seqs = read_gfa(para, DAG, para->inc_fp.c_str());

  int exist_seq_num = seqs.size();
  try {
    if (!path.empty()) readFile(seqs, path.c_str());
    // std::cerr << exist_seq_num << " " << seqs.size() << " " << DAG->node.size() << "\n";
  }
  catch (const std::exception& e) {
    std::cerr << "error read file: " << e.what() << std::endl;
    return 1;
  }

  if (para->verbose && para->enable_seeding) std::cerr << "collect minimizer" << "\n";
  minimizer_t* mm = new minimizer_t(para, seqs);
  if (para->verbose && para->progressive_poa) std::cerr << "build guide tree" << "\n";
  if (para->inc_fp.empty() && para->progressive_poa) {
    if (para->verbose) std::cerr << "progressive" << "\n";
    mm->get_guide_tree(para);
  }
  const std::vector<int>& ord = mm->ord;

  int rid;
  if (para->verbose) std::cerr << "poa" << "\n";
  aligned_buff_t* mpool = new aligned_buff_t[para->thread];
  for (int i = exist_seq_num; i < seqs.size(); i++) {  //seqs.size()
    rid = ord[i];
    if (para->verbose && i % 10 == 0) {
      std::cerr << "[" << i << "/" << seqs.size() << "]" << "\n";
    }
    if (para->verbose) std::cerr << "aligment" << "\n";
    std::vector<res_t> res = alignment(para, DAG, mm, rid, seqs[rid].seq, mpool);
    if (para->verbose) std::cerr << "add path" << "\n";
    DAG->add_path(para->m, i, res, 1);
    if (para->verbose) std::cerr << "topsort" << "\n";
    DAG->topsort(para, i + 1 == seqs.size());
  }

  if (para->verbose) std::cerr << "out_put" << "\n";
  if (para->result == 0) DAG->output_consensus();
  else if (para->result == 1) DAG->output_rc_msa(para, mm->rid_to_ord, seqs);
  else if (para->result == 2) DAG->output_gfa(mm->rid_to_ord, seqs);

  // delete
  delete[] mpool;
  mpool = nullptr;  // 防止后续误用
  delete para;
  para = nullptr;  // 防止后续误用
  delete DAG;
  DAG = nullptr;  // 防止后续误用
  delete mm;
  mm = nullptr;  // 防止后续误用
  return 0;
}