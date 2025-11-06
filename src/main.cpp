#include <iostream>
#include "cxxopts.hpp"
#include "handle_input.h"
#include "parameter.h"
#include "graph.h"
#include "align.h"
#include <immintrin.h>
#include "mem_alloc_utils.h"
#include "ThreadPool.h"
#include <algorithm>
#include <numeric>
#include "minimizer.h"
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
    ("t,thread", "thread number", cxxopts::value<int>()->default_value("0"))
    ("b,band_b", "band arg", cxxopts::value<int>()->default_value("100"))
    ("f,band_f", "band arg", cxxopts::value<int>()->default_value("40"))
    ("B,ab_band", "adpative band arg", cxxopts::value<bool>()->default_value("false"))
    ("S,seeding", "enable minimizer-based seeding and anchoring", cxxopts::value<bool>()->default_value("false"))
    ("W,poa_w", "the minimum distance between adjacent anchors", cxxopts::value<int>()->default_value("500"))
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
    para->match = result["match"].as<int>();
    para->mismatch = result["mismatch"].as<int>();
    para->gap_open1 = result["gap_open"].as<int>();
    para->gap_ext1 = result["gap_ext"].as<int>();
    para->b = result["band_b"].as<int>();
    para->f = result["band_f"].as<int>();
    para->ab_band = result["ab_band"].as<bool>();
    para->k = result["k_mer"].as<int>();
    para->mm_w = result["mm_w"].as<int>();
    para->progressive_poa = result["progressive_poa"].as<bool>();
    para->enable_seeding = result["seeding"].as<bool>();
    thread = result["thread"].as<int>();
    sample_num = result["sample_num"].as<int>();
    if (sample_num < 0) sample_num * -1;
    para->result = result["result"].as<int>();
    para->verbose = result["verbose"].as<int>();
    para->poa_w = result["poa_w"].as<int>();
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
  // std::cerr << DAG->node.size() << "\n";
  // return 0;
  int exist_seq_num = seqs.size();
  try {
    if (!path.empty()) readFile(seqs, path.c_str());
    // std::cerr << exist_seq_num << " " << seqs.size() << " " << DAG->node.size() << "\n";
  }
  catch (const std::exception& e) {
    std::cerr << "error read file: " << e.what() << std::endl;
    return 1;
  }
  // std::cerr << seqs.size() << "\n";
  // seqs.resize(660);
  // handle alignment 
  // std::sort(ord.begin(), ord.end(), [&](const int& i, const int& j) {
  //   return  seqs[i].seq.size() > seqs[j].seq.size();
  // });
  if (para->verbose && para->enable_seeding) std::cerr << "collect minimizer" << "\n";
  minimizer_t* mm = new minimizer_t(para, seqs);
  if (para->verbose && para->progressive_poa) std::cerr << "build guide tree" << "\n";
  if (para->inc_fp.empty() && para->progressive_poa) {
    std::cerr << "progressive" << "\n";
    mm->get_guide_tree(para);
    // std::reverse(ord.begin(), ord.end());
    // for (int i = 0; i < seqs.size(); i++) {
    //   std::cerr << ord[i] << "\n";
    // }
  }
  const std::vector<int>& ord = mm->ord;
  // std::cerr << mm.mm_v.n << "\n";

  int rid;
  // for (int i = 0; i < seqs.size(); i++) {
  //   std::cout << i << " " << ord[i] << " " << seqs[ord[i]].seq.size() << "\n";
  // }
  if (para->verbose) std::cerr << "poa" << "\n";
  // seqs[0].seq = "TTGCCCTT";
  // seqs[1].seq = "CCAATTTT";
  // seqs[2].seq = "TGCT";
  if (!thread) { // 
    aligned_buff_t mpool;
    for (int i = exist_seq_num; i < seqs.size(); i++) {  //seqs.size()
      rid = ord[i];
      if (para->verbose && i % 10 == 0) {
        std::cerr << "[" << i << "/" << seqs.size() << "]" << "\n";
      }

      // std::string tseq;
      // tseq += char26_table['N'];
      // for (int j = 0; j < seqs[rid].seq.size(); j++) {
      //   tseq += char26_table[seqs[rid].seq[j]];
      // }
      // std::cerr << "poa" << "\n";
      // POA_SIMD(para, DAG, tseq);
      // std::vector<res_t> res = POA(para, DAG, tseq);
      // std::vector<res_t> res = POA_SIMD(para, DAG, tseq);
      // std::vector<res_t> res = abPOA(para, DAG, mm, rid, tseq, &mpool);
      if (para->verbose) std::cerr << "aligment" << "\n";
      std::vector<res_t> res = alignment(para, DAG, mm, rid, seqs[rid].seq, &mpool);
      // return 0;
      // std::cerr << "add_path" << "\n";
      if (para->verbose) std::cerr << "add path" << "\n";
      DAG->add_path(para->m, i, res);
      // std::cerr << "topsort" << "\n";
      if (para->verbose) std::cerr << "topsort" << "\n";
      DAG->topsort(para, i + 1 == seqs.size());
      // std::cout << i << " " << DAG->rank.size() << "\n";
    }
    // handle output 
    // DAG->output_rc_msa(mm->rid_to_ord, seqs);
    // std::cerr << minl << " " << maxl << '\n';
  }
  else {
    ThreadPool pool(thread);
    aligned_buff_t* mpool = new aligned_buff_t[thread];
    // std::cerr << "thread:" << thread << "\n";
    std::vector<std::future<std::vector<res_t>> > results;
    // sample_num = 0;
    for (int i = exist_seq_num; i < seqs.size() && i - exist_seq_num < sample_num; i++) { // ensure parallel before DAG have the enough sample seq
      rid = ord[i];
      if (para->verbose && i % 10 == 0) {
        std::cerr << "[" << i << "/" << seqs.size() << "]" << "\n";
      }
      // std::string tseq;
      // tseq += char26_table['N'];
      // for (int j = 0; j < seqs[rid].seq.size(); j++) {
      //   tseq += char26_table[seqs[rid].seq[j]];
      // }
      // std::cerr << "poa" << "\n";
      // POA_SIMD(para, DAG, tseq);
      // std::vector<res_t> res = POA(para, DAG, tseq);
      // std::vector<res_t> res = POA_SIMD(para, DAG, tseq);
      // std::vector<res_t> res = abPOA(para, DAG, mm, rid, tseq, &mpool[0]);
      std::vector<res_t> res = alignment(para, DAG, mm, rid, seqs[rid].seq, &mpool[0]);
      // return 0;
      // std::cerr << "add_path" << "\n";
      DAG->add_path(para->m, i, res);
      // std::cerr << "topsort" << "\n";
      DAG->topsort(para, i + 1 == seqs.size());
    }
    for (int i = exist_seq_num + sample_num; i < seqs.size(); i += thread) {  //seqs.size()
      results.clear();
      for (int j = 0; j < thread && i + j < seqs.size(); j++) {
        rid = ord[i + j];
        // const char* seq_i = seqs[i].seq.c_str();
        const std::string& seq_i = seqs[rid].seq;
        aligned_buff_t* cur_mpool = &mpool[j];
        results.emplace_back(
          pool.enqueue([para, DAG, mm, rid, &seq_i, cur_mpool] {
          // std::string tseq;
          // tseq += char26_table['N'];
          // for (int j = 0; j < seq_i.size(); j++) {
          //   tseq += char26_table[seq_i[j]];
          // }
          // std::vector<res_t> res = abPOA(para, DAG, mm, rid, tseq, cur_mpool);
          std::vector<res_t> res = alignment(para, DAG, mm, rid, seq_i, cur_mpool);
          return res;
        }));
      }
      std::vector<std::vector<res_t>> res;
      if (para->verbose) {
        std::cerr << "[" << i << "/" << seqs.size() << "]" << "\n";
      }
      for (auto&& result : results) {
        res.emplace_back(result.get());
      }

      // std::cerr << "add_path" << "\n";
      int node_num = DAG->node.size();
      for (int j = 0; j < res.size(); j++) {
        rid = ord[i + j];
        DAG->add_path(para->m, i + j, res[j], node_num);
      }
      // std::cerr << "topsort:" << "\n";
      DAG->topsort(para, i + thread >= seqs.size());
      // std::cerr << "poa" << "\n";
      // POA_SIMD(para, DAG, tseq);
      // std::vector<res_t> res = POA(para, DAG, tseq);
      // std::vector<res_t> res = POA_SIMD(para, DAG, tseq);
      // std::cerr << "topsort" << "\n";
      // std::cout << i << " " << DAG->rank.size() << "\n";
    }
    // handle output 
    // std::cerr << "output" << '\n';
    // std::cerr << "delete" << "\n";
    delete[] mpool;
    // std::cerr << "finish" << "\n";
  }
  // std::cerr << "out_put" << "\n";
  if (para->result == 0) DAG->output_consensus();
  else if (para->result == 1) DAG->output_rc_msa(para, mm->rid_to_ord, seqs);
  else if (para->result == 2) DAG->output_gfa(mm->rid_to_ord, seqs);


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
  delete DAG;
  DAG = nullptr;  // 防止后续误用
  delete mm;
  mm = nullptr;  // 防止后续误用
  return 0;
}