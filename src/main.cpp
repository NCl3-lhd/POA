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
    ("i,input", "input path", cxxopts::value<std::string>())
    ("m,mat_fp", "match file path", cxxopts::value<std::string>())
    ("M,match", "match sorce", cxxopts::value<int>()->default_value("2"))
    ("X,mismatch", "mismatch sorce", cxxopts::value<int>()->default_value("-4"))
    ("O,gap_open", "gap_open sorce", cxxopts::value<int>()->default_value("-4"))
    ("E,gap_ext", "gap_ext sorce", cxxopts::value<int>()->default_value("-2"))
    ("t,thread", "thread number", cxxopts::value<int>()->default_value("0"))
    ("b,band_b", "band arg", cxxopts::value<int>()->default_value("50"))
    ("f,band_f", "band arg", cxxopts::value<int>()->default_value("100"))
    ("k,k_mer", "k_mer lenth", cxxopts::value<int>()->default_value("19")) //19
    ("w,window", "k_mer_window lenth", cxxopts::value<int>()->default_value("10"))
    ("s,sample_num", "sample_num", cxxopts::value<int>()->default_value("50"))
    ("p,progressive_poa", "is progressive_poa", cxxopts::value<bool>()->default_value("false"))
    ("h,help", "Print usage")
    ;
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
    para->match = result["match"].as<int>();
    para->mismatch = result["mismatch"].as<int>();
    para->gap_open1 = result["gap_open"].as<int>();
    para->gap_ext1 = result["gap_ext"].as<int>();
    para->b = result["band_b"].as<int>();
    para->f = result["band_f"].as<int>();
    para->k = result["k_mer"].as<int>();
    para->w = result["window"].as<int>();
    para->progressive_poa = result["progressive_poa"].as<bool>();
    thread = result["thread"].as<int>();
    sample_num = result["sample_num"].as<int>();
    if (sample_num < 0) sample_num * -1;
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
  // std::cerr << seqs.size() << "\n";
  // seqs.resize(53);
  // handle alignment 
  graph* DAG = new graph();
  std::vector<int> ord(seqs.size());
  std::iota(ord.begin(), ord.end(), 0);
  // std::sort(ord.begin(), ord.end(), [&](const int& i, const int& j) {
  //   return  seqs[i].seq.size() > seqs[j].seq.size();
  // });
  minimizer_t* mm = nullptr;
  if (para->progressive_poa) {
    std::cerr << "progressive_poa" << "\n";
    mm = new minimizer_t(para, seqs);
    ord = mm->get_guide_tree(para);
    // std::reverse(ord.begin(), ord.end());
    // for (int i = 0; i < seqs.size(); i++) {
    //   std::cerr << ord[i] << "\n";
    // }
  }
  // std::cerr << mm.mm_v.n << "\n";


  // for (int i = 0; i < seqs.size(); i++) {
  //   std::cout << i << " " << ord[i] << " " << seqs[ord[i]].seq.size() << "\n";
  // }
  int seq_id = ord[0];
  DAG->init(para, seq_id, seqs[seq_id].seq);
  // seqs[0].seq = "TTGCCCTT";
  // seqs[1].seq = "CCAATTTT";
  // seqs[2].seq = "TGCT";
  if (!thread) { // 
    aligned_buff_t mpool;
    for (int i = 1; i < seqs.size(); i++) {  //seqs.size()
      seq_id = ord[i];
      if (i % 10 == 0) {
        std::cerr << "[" << i << "/" << seqs.size() << "]" << "\n";
      }
      std::string tseq;
      tseq += char26_table['N'];
      for (int j = 0; j < seqs[seq_id].seq.size(); j++) {
        tseq += char26_table[seqs[seq_id].seq[j]];
      }
      // std::cerr << "poa" << "\n";
      // POA_SIMD(para, DAG, tseq);
      // std::vector<res_t> res = POA(para, DAG, tseq);
      // std::vector<res_t> res = POA_SIMD(para, DAG, tseq);
      std::vector<res_t> res = para->f ? abPOA(para, DAG, mm, tseq, &mpool) : POA_SIMD_ORIGIN(para, DAG, tseq, &mpool);
      // return 0;
      // std::cerr << "add_path" << "\n";
      DAG->add_path(para->m, seq_id, res);
      // std::cerr << "topsort" << "\n";
      DAG->topsort(i + 1 == seqs.size(), para->f);
      // std::cout << i << " " << DAG->rank.size() << "\n";
    }
    // handle output 
    DAG->output_rc_msa(seqs);
    // std::cerr << minl << " " << maxl << '\n';
  }
  else {
    ThreadPool pool(thread);
    aligned_buff_t* mpool = new aligned_buff_t[thread];
    // std::cerr << "thread:" << thread << "\n";
    std::vector<std::future<std::vector<res_t>> > results;
    // sample_num = 0;
    for (int i = 1; i < seqs.size() && i <= sample_num; i++) { // ensure parallel before DAG have the enough sample seq
      seq_id = ord[i];
      if (i % 10 == 0) {
        std::cerr << "[" << i << "/" << seqs.size() << "]" << "\n";
      }
      std::string tseq;
      tseq += char26_table['N'];
      for (int j = 0; j < seqs[seq_id].seq.size(); j++) {
        tseq += char26_table[seqs[seq_id].seq[j]];
      }
      // std::cerr << "poa" << "\n";
      // POA_SIMD(para, DAG, tseq);
      // std::vector<res_t> res = POA(para, DAG, tseq);
      // std::vector<res_t> res = POA_SIMD(para, DAG, tseq);
      std::vector<res_t> res = para->f ? abPOA(para, DAG, mm, tseq, &mpool[0]) : POA_SIMD_ORIGIN(para, DAG, tseq, &mpool[0]);
      // return 0;
      // std::cerr << "add_path" << "\n";
      DAG->add_path(para->m, seq_id, res);
      // std::cerr << "topsort" << "\n";
      DAG->topsort(i + 1 == seqs.size(), para->f);
    }
    for (int i = sample_num + 1; i < seqs.size(); i += thread) {  //seqs.size()
      results.clear();
      for (int j = 0; j < thread && i + j < seqs.size(); j++) {
        seq_id = ord[i + j];
        if (j % 3 == 0) {
          std::cerr << "[" << i + j << "/" << seqs.size() << "]" << "\n";
        }
        // const char* seq_i = seqs[i].seq.c_str();
        const std::string& seq_i = seqs[seq_id].seq;
        aligned_buff_t* cur_mpool = &mpool[j];
        results.emplace_back(
          pool.enqueue([para, DAG, mm, &seq_i, cur_mpool] {
          std::string tseq;
          tseq += char26_table['N'];
          for (int j = 0; j < seq_i.size(); j++) {
            tseq += char26_table[seq_i[j]];
          }
          std::vector<res_t> res = para->f ? abPOA(para, DAG, mm, tseq, cur_mpool) : POA_SIMD_ORIGIN(para, DAG, tseq, cur_mpool);
          return res;
        }));
      }
      // std::cerr << "add_path" << "\n";
      std::vector<std::vector<res_t>> res;
      for (auto&& result : results) {
        res.emplace_back(result.get());
      }
      int node_num = DAG->node.size();
      for (int j = 0; j < res.size(); j++) {
        seq_id = ord[i + j];
        DAG->add_path(para->m, seq_id, res[j], node_num);
      }
      // std::cerr << "topsort:" << "\n";
      DAG->topsort(i + thread >= seqs.size(), para->f);

      // std::cerr << "poa" << "\n";
      // POA_SIMD(para, DAG, tseq);
      // std::vector<res_t> res = POA(para, DAG, tseq);
      // std::vector<res_t> res = POA_SIMD(para, DAG, tseq);
      // std::cerr << "topsort" << "\n";
      // std::cout << i << " " << DAG->rank.size() << "\n";
    }
    // handle output 
    // std::cerr << "output" << '\n';
    DAG->output_rc_msa(seqs);
    // std::cerr << "delete" << "\n";
    delete[] mpool;
    // std::cerr << "finish" << "\n";
  }



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