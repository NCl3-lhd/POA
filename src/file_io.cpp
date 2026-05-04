#include <zlib.h>
#include <stdexcept> 
#include <stdint.h>
#include <unistd.h>
#include "file_io.h"
#include "kseq.h"
#include "sequence.h" 
#include "utils.h"
#include <stdlib.h>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <cctype>
#include <string>

constexpr int INF = 0x3f3f3f3f; // 0x3f3f3f3f

KSEQ_INIT(gzFile, gzread)

// handle input file
void readFile(para_t* para, std::vector<seq_t>& seqs, const char* path) {
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
  // std::vector<seq_t> seqs;
  while ((l = kseq_read(seq)) >= 0) {
    // printf("name: %s\n", seq->name.s); // 打印序列名称
    seq_t tseq;
    tseq.name = seq->name.s;
    if (seq->comment.l) tseq.comment = seq->comment.s;// 序列注释
    tseq.seq = seq->seq.s;
    for (int i = 0; i < tseq.seq.size(); i++) {
      para->isRNA |= (tseq.seq[i] | 32) == 'u';
    }
    if (seq->qual.l)  tseq.qual = seq->qual.s;  // 质量分数
    seqs.emplace_back(tseq);
  }
  // 释放 kseq 结构体
  kseq_destroy(seq);
  // 关闭文件
  gzclose(fp);
  return;
}
// handle input arg
void initPara(para_t* para) {
  para->mm_filter_ratio = 0.25f;
  //handle incorrect para
  if (para->match < 0) para->m *= -1;
  if (para->mismatch > 0) para->mismatch *= -1;
  if (para->gap_open1 > 0) para->gap_open1 *= -1;
  if (para->gap_ext1 > 0) para->gap_ext1 *= -1;
  if (para->mat_fp.empty()) {
    para->m = 5; // default m = Nucleotide num
    int m = para->m;
    para->mat.resize(m * m, 0); //like HOXD70.mtx
    for (int i = 0; i < m - 1; i++) {
      for (int j = 0; j < m - 1; j++) {
        para->mat[i * m + j] = i == j ? para->match : para->mismatch;
      }
    }
    para->mat[(m - 1) * m + m - 1] = para->match;
  }
  else {
    std::ifstream mtx_file(para->mat_fp);
    if (!mtx_file.is_open()) {
      throw std::runtime_error("Error opening matrix file " + para->mat_fp);
    }
    std::string line;
    while (std::getline(mtx_file, line) && (!line.empty() && line[0] == '#')) {}
    if (line.empty()) {
      throw std::runtime_error("Matrix file " + para->mat_fp + " has no header line");
    }
    std::istringstream header(line);
    std::vector<char> alphabet;
    std::string token;
    while (header >> token) {
      if (token.size() == 1) alphabet.push_back(token[0]);
    }
    int m = alphabet.size();
    if (m == 0) {
      throw std::runtime_error("Matrix file " + para->mat_fp + " header is empty");
    }
    para->m = m;
    std::vector<int> col_to_idx(m);
    for (int i = 0; i < m; i++) {
      unsigned char c = alphabet[i];
      col_to_idx[i] = m > 5 ? aa26_table[c] : nt4_table[c];
    }
    para->mat.resize(m * m, 0);
    while (std::getline(mtx_file, line)) {
      if (line.empty() || line[0] == '#') continue;
      std::istringstream data(line);
      std::string row_label;
      if (!(data >> row_label)) continue;
      if (row_label.size() != 1) {
        throw std::runtime_error("Invalid row label in matrix file " + para->mat_fp + ": " + row_label);
      }
      unsigned char row_char = row_label[0];
      int row_idx = m > 5 ? aa26_table[row_char] : nt4_table[row_char];
      int score;
      for (int j = 0; j < m; j++) {
        if (!(data >> score)) {
          throw std::runtime_error("Not enough columns in matrix file " + para->mat_fp + " row " + row_label);
        }
        para->mat[row_idx * m + col_to_idx[j]] = score;
      }
    }
    mtx_file.close();
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
int gfa_parse_S(para_t* para, graph* DAG, seg_seq_t* segs, char* s) {
  if (s[1] != '\t' || s[2] == '\0') return -1;
  char* end_s, * start_s, * seg = 0;
  int i, seg_len, seg_name_len, is_ok = 0;
  char* seg_name = 0;

  for (i = 0, end_s = start_s = s + 2;; ++end_s) {
    if (*end_s == '\0' || *end_s == '\t') {
      int c = *end_s;
      *end_s = '\0';
      if (i == 0) {
        seg_name = start_s;
        seg_name_len = end_s - start_s;
      }
      else if (i == 1) {
        seg = start_s;
        seg_len = end_s - start_s;
        is_ok = 1;
        break;
      }
      if (c == '\0') break;
      ++i, start_s = end_s + 1;
    }
  }

  if (is_ok) {
    seg_seq_realloc(segs);
    kputsn(seg_name, seg_name_len, segs->name + segs->n);
    kputsn(seg, seg_len, segs->seq + segs->n);
    int absent, id;
    khint_t pos = kh_put(str, segs->h, segs->name[segs->n].s, &absent);
    if (absent) {
      kh_val(segs->h, pos) = segs->n;
      // add seg_seq into DAG
      if (absent) { // add node for seg_seq
        kv_resize(uint32_t, NULL, segs->start_id, segs->start_id.n + seg_len);
        // kv_push(uint32_t, NULL, segs->start_id, id);
        for (i = 0; i < seg_len; ++i) {
          id = DAG->add_node(para, seg[i]);
          // id = abpoa_add_graph_node(abg, ab_char26_table[(int)(seg->s[i])]);
          if (i == 0) kv_push(uint32_t, NULL, segs->start_id, id);
          if (i == seg_len - 1) kv_push(uint32_t, NULL, segs->end_id, id);
        }
      }
      //
    }
    else err_fatal(__func__, "Duplicated chromosome: \"%s\".", seg_name);
    ++segs->n;
  }
  else err_fatal(__func__, "Error: no seq in GFA segment line (%s).", seg_name);
  return 0;

}
int gfa_parse_P(para_t *para, graph *DAG, std::vector<seq_t> &seqs, seg_seq_t *segs, char *s, PathWriter *writer) {
  if (s[1] != '\t' || s[2] == '\0') return -1;
  char* end_s, * start_s, * path = 0;
  int i, is_ok = 0, is_rc = -1;
  char* path_name = 0; int path_name_len = 0;
  kstring_t* seg_seq, * seg_name;

  for (i = 0, end_s = start_s = s + 2;; ++end_s) {
    if (*end_s == 0 || *end_s == '\t') {
      int c = *end_s;
      *end_s = 0;
      if (i == 0) {
        path_name = start_s;
        path_name_len = end_s - start_s;
      }
      else if (i == 1) {
        path = start_s;
        is_ok = 1;
        break;
      }
      if (c == 0) break;
      ++i, start_s = end_s + 1;
    }
  }

  if (is_ok) {
    std::vector<int> path_node_ids;
    char* deli_s, * info_s, * _seg_name; khint_t pos, seg_pos; int absent;
    int id, in_id = -1, out_id = -1, last_id = 0, next_id = 1;
    int curPos = 0;
    for (deli_s = info_s = path; ; ++deli_s) {
      if (*deli_s == '+') {
        if (is_rc == 1) err_fatal(__func__, "Error: path has both \'+\' and \'-\' seg. (%s)", path_name);
        is_rc = 0; *deli_s = 0; _seg_name = info_s;
        seg_pos = kh_get(str, segs->h, _seg_name);
        if (seg_pos == kh_end(segs->h)) err_fatal(__func__, "Error: seg (%s) not exist.", info_s);
        seg_name = segs->name + kh_val(segs->h, seg_pos);
        seg_seq = segs->seq + kh_val(segs->h, seg_pos);

        // check if seg already exist
        // pos = kh_put(abstr, seg_name2in_id, seg_name->s, &absent);
        in_id = segs->start_id.a[kh_val(segs->h, seg_pos)];
        out_id = segs->end_id.a[kh_val(segs->h, seg_pos)];
        // if (absent) { // add node for seg_seq
        //   for (i = 0; i < (int)seg_seq->l; ++i) {
        //     id = abpoa_add_graph_node(abg, ab_char26_table[(int)(seg_seq->s[i])]);
        //     if (i == 0) in_id = id;
        //     if (i == (int)seg_seq->l - 1) out_id = id;
        //   }
        //   kh_val(seg_name2in_id, pos) = in_id;
        //   pos = kh_put(str, seg_name2out_id, seg_name->s, &absent);
        //   kh_val(seg_name2out_id, pos) = out_id;
        // }
        // else {
        //   in_id = kh_val(seg_name2in_id, pos);
        //   pos = kh_put(str, seg_name2out_id, seg_name->s, &absent);
        //   out_id = kh_val(seg_name2out_id, pos);
        // }
        // add edge
        // std::cerr << last_id << " " << in_id << "\n";
        // abpoa_add_graph_edge(abg, last_id, in_id, 1, 1, add_read_id, 0, p_i, read_ids_n, p_n);
        DAG->add_adj(last_id, in_id);
        if (in_id < out_id) {
          for (i = 0; i < out_id - in_id; ++i) {
            // std::cerr << in_id + i << " " << in_id + i + 1 << "\n";
            DAG->add_adj(in_id + i, in_id + i + 1);
          }
          for (i = in_id; i <= out_id; i++) {
            path_node_ids.emplace_back(i);
          }
        }
        else if (in_id > out_id) err_fatal(__func__, "Error: in_id (%d) > out_id (%d).", in_id, out_id);

        last_id = out_id;
        info_s = deli_s + 2;
      }
      else if (*deli_s == '-') {
        if (is_rc == 0) err_fatal(__func__, "Error: path has both \'+\' and \'-\' seg. (%s)", path_name);
        is_rc = 1; *deli_s = 0; _seg_name = info_s;
        seg_pos = kh_get(str, segs->h, _seg_name);
        if (seg_pos == kh_end(segs->h)) err_fatal(__func__, "Error: seg (%s) not exist.", info_s);
        seg_name = segs->name + kh_val(segs->h, seg_pos);
        seg_seq = segs->seq + kh_val(segs->h, seg_pos);

        // check if seg exist
        in_id = segs->start_id.a[kh_val(segs->h, seg_pos)];
        out_id = segs->end_id.a[kh_val(segs->h, seg_pos)];
        // pos = kh_put(str, seg_name2in_id, seg_name->s, &absent);
        // if (absent) { // add node for seg_seq
        //   for (i = 0; i < (int)seg_seq->l; ++i) {
        //     id = abpoa_add_graph_node(abg, ab_char26_table[(int)(seg_seq->s[i])]);
        //     if (i == 0) in_id = id;
        //     if (i == (int)seg_seq->l - 1) out_id = id;
        //   }
        //   kh_val(seg_name2in_id, pos) = in_id;
        //   pos = kh_put(str, seg_name2out_id, seg_name->s, &absent);
        //   kh_val(seg_name2out_id, pos) = out_id;
        // }
        // else {
        //   in_id = kh_val(seg_name2in_id, pos); out_id = kh_val(seg_name2out_id, pos);
        // }

        // add edge

        DAG->add_adj(out_id, next_id);

        // abpoa_add_graph_edge(abg, out_id, next_id, 1, 1, add_read_id, 0, p_i, read_ids_n, p_n);
        if (in_id < out_id) {
          for (i = 0; i < out_id - in_id; ++i) {
            DAG->add_adj(in_id + i, in_id + i + 1);
          }
          for (i = out_id; i >= in_id; --i) {
            path_node_ids.emplace_back(i);
          }
          // abpoa_add_graph_edge(abg, in_id + i, in_id + i + 1, 1, 1, add_read_id, 0, p_i, read_ids_n, p_n);
        }
        else if (in_id > out_id) err_fatal(__func__, "Error: in_id (%d) > out_id (%d).", in_id, out_id);

        next_id = in_id;
        info_s = deli_s + 2;
      }
      else if (*deli_s == 0 || *deli_s == '\t') break;
    }
    if (is_rc) {
      // abpoa_add_graph_edge(abg, ABPOA_SRC_NODE_ID, next_id, 1, 1, add_read_id, 0, p_i, read_ids_n, p_n);
      DAG->add_adj(0, next_id);
    }
    else {
      // abpoa_add_graph_edge(abg, last_id, ABPOA_SINK_NODE_ID, 1, 1, add_read_id, 0, p_i, read_ids_n, p_n);
      DAG->add_adj(last_id, 1);
    }
    // set abs
    // abpoa_cpy_str(abs->name + abs->n_seq, walk_name, walk_name_len);
    seq_t tseq;
    tseq.name = path_name;
    if (para->result) writer->write_path(seqs.size(), path_node_ids);
    seqs.emplace_back(tseq);
    // abs->is_rc[abs->n_seq] = is_rc; abs->n_seq++;
  }
  else err_fatal(__func__, "Error: no path in GFA path line (%s).", path_name);
  return 0;
}
// int gfa_parse_W(para_t* para, graph* DAG, std::vector<seq_t>& seqs, seg_seq_t* segs, char* s) {
//   if (s[1] != '\t' || s[2] == '\0') return -1;
//   char* end_s, * start_s, * path = 0;
//   int i, is_ok = 0, is_rc = -1;
//   char* walk_name = 0; int walk_name_len = 0;
//   kstring_t* seg_seq, * seg_name;
//   int start_pos, end_pos;
//   // std::cerr << s << "\n";
//   for (i = 0, end_s = start_s = s + 2;; ++end_s) {
//     if (*end_s == '\0' || *end_s == '\t') {
//       int c = *end_s;
//       *end_s = '\0';
//       if (i == 0) {
//         // walk id
//       }
//       else if (i == 1) {

//       }
//       else if (i == 2) {
//         walk_name = start_s;
//         walk_name_len = end_s - start_s;
//       }
//       else if (i == 3) {
//         start_pos = atoi(start_s);
//       }
//       else if (i == 4) {
//         end_pos = atoi(start_s);
//       }
//       else if (i == 5) {
//         path = start_s;
//         is_ok = 1;
//         break;
//       }
//       if (c == 0) break;
//       ++i, start_s = end_s + 1;
//     }
//   }
//   // std::cerr << walk_name << " " << path << "\n";
//   // std::cerr << "ok" << "\n";
//   // std::cerr << start_pos << " " << end_pos << "\n";
//   std::vector<int> path_node;
//   if (is_ok) {
//     char* end_s, * start_s, * _seg_name; khint_t pos, seg_pos; int absent;
//     int id, in_id = -1, out_id = -1, last_id = 0, next_id = 1;
//     int curPos = 0;
//     for (end_s = start_s = path + 1; ; ++end_s) {
//       if (*end_s == '>') {
//         if (is_rc == 1) err_fatal(__func__, "Error: walk has both \'>\' and \'<\' seg. (%s)", walk_name);
//         is_rc = 0; *end_s = '\0'; _seg_name = start_s;
//         // std::cerr << _seg_name << " ";
//         // if (_seg_name == "s211") {
//         //   exit(1);
//         // }
//         seg_pos = kh_get(str, segs->h, _seg_name);
//         if (seg_pos == kh_end(segs->h)) err_fatal(__func__, "Error: seg (%s) not exist.", start_s);
//         seg_name = segs->name + kh_val(segs->h, seg_pos);
//         seg_seq = segs->seq + kh_val(segs->h, seg_pos);

//         // check if seg already exist
//         in_id = segs->start_id.a[kh_val(segs->h, seg_pos)];
//         out_id = segs->end_id.a[kh_val(segs->h, seg_pos)];
//         // pos = kh_put(str, seg_name2in_id, seg_name->s, &absent);
//         // if (absent) { // add node for seg_seq
//         //   for (i = 0; i < (int)seg_seq->l; ++i) {
//         //     id = abpoa_add_graph_node(abg, ab_char26_table[(int)(seg_seq->s[i])]);
//         //     if (i == 0) in_id = id;
//         //     if (i == (int)seg_seq->l - 1) out_id = id;
//         //   }
//         //   kh_val(seg_name2in_id, pos) = in_id;
//         //   pos = kh_put(str, seg_name2out_id, seg_name->s, &absent);
//         //   kh_val(seg_name2out_id, pos) = out_id;
//         // }
//         // else {
//         //   in_id = kh_val(seg_name2in_id, pos);
//         //   pos = kh_put(str, seg_name2out_id, seg_name->s, &absent);
//         //   out_id = kh_val(seg_name2out_id, pos);
//         // }
//         // add edge
//         // std::cerr << in_id << " " << out_id << "\n";
//         // DAG->add_adj(seqs.size(), last_id, in_id, curPos++);
//         path_node.push_back(in_id);
//         // std::cerr << last_id << " " << in_id << "\n";
//         // abpoa_add_graph_edge(abg, last_id, in_id, 1, 1, add_read_id, 0, p_i, read_ids_n, p_n);
//         if (in_id < out_id) {
//           for (i = 0; i < out_id - in_id; ++i) {
//             // std::cerr << in_id + i << " " << in_id + i + 1 << "\n";
//             path_node.push_back(in_id + i + 1);
//             // DAG->add_adj(seqs.size(), in_id + i, in_id + i + 1, curPos++);
//           }
//         }
//         else if (in_id > out_id) err_fatal(__func__, "Error: in_id (%d) > out_id (%d).", in_id, out_id);

//         last_id = out_id;
//         start_s = end_s + 1;
//       }
//       else if (*end_s == '<') { // is not correct, not support <
//         if (is_rc == 0) err_fatal(__func__, "Error: walk has both \'>\' and \'<\' seg. (%s)", walk_name);
//         is_rc = 1; *end_s = '\0'; _seg_name = start_s;
//         seg_pos = kh_get(str, segs->h, _seg_name);
//         if (seg_pos == kh_end(segs->h)) err_fatal(__func__, "Error: seg (%s) not exist.", start_s);
//         seg_name = segs->name + kh_val(segs->h, seg_pos);
//         seg_seq = segs->seq + kh_val(segs->h, seg_pos);

//         // check if seg exist
//         in_id = segs->start_id.a[seg_pos];
//         out_id = segs->end_id.a[seg_pos];
//         // pos = kh_put(str, seg_name2in_id, seg_name->s, &absent);
//         // if (absent) { // add node for seg_seq
//         //   for (i = 0; i < (int)seg_seq->l; ++i) {
//         //     id = abpoa_add_graph_node(abg, ab_char26_table[(int)(seg_seq->s[i])]);
//         //     if (i == 0) in_id = id;
//         //     if (i == (int)seg_seq->l - 1) out_id = id;
//         //   }
//         //   kh_val(seg_name2in_id, pos) = in_id;
//         //   pos = kh_put(str, seg_name2out_id, seg_name->s, &absent);
//         //   kh_val(seg_name2out_id, pos) = out_id;
//         // }
//         // else {
//         //   in_id = kh_val(seg_name2in_id, pos); out_id = kh_val(seg_name2out_id, pos);
//         // }

//         // add edge
//         DAG->add_adj(seqs.size(), out_id, next_id, curPos++);
//         // abpoa_add_graph_edge(abg, out_id, next_id, 1, 1, add_read_id, 0, p_i, read_ids_n, p_n);
//         if (in_id < out_id) {
//           for (i = 0; i < out_id - in_id; ++i)
//             DAG->add_adj(seqs.size(), in_id + i, in_id + i + 1, curPos++);
//           // abpoa_add_graph_edge(abg, in_id + i, in_id + i + 1, 1, 1, add_read_id, 0, p_i, read_ids_n, p_n);
//         }
//         else if (in_id > out_id) err_fatal(__func__, "Error: in_id (%d) > out_id (%d).", in_id, out_id);

//         next_id = in_id;
//         start_s = end_s + 1;
//       }
//       else if (*end_s == '\0' || *end_s == '\t') break;
//     }

//     if (is_rc) {
//       // abpoa_add_graph_edge(abg, ABPOA_SRC_NODE_ID, next_id, 1, 1, add_read_id, 0, p_i, read_ids_n, p_n);
//       DAG->add_adj(seqs.size(), 0, next_id, curPos++);
//     }
//     else {
//       //abpoa_add_graph_edge(abg, last_id, ABPOA_SINK_NODE_ID, 1, 1, add_read_id, 0, p_i, read_ids_n, p_n);
//       // std::cerr << start_s << " ";
//       seg_pos = kh_get(str, segs->h, start_s);
//       if (seg_pos == kh_end(segs->h)) err_fatal(__func__, "Error: seg (%s) not exist.", start_s);
//       seg_seq = segs->seq + kh_val(segs->h, seg_pos);

//       // check if seg already exist
//       in_id = segs->start_id.a[kh_val(segs->h, seg_pos)];
//       out_id = segs->end_id.a[kh_val(segs->h, seg_pos)];
//       // DAG->add_adj(seqs.size(), last_id, in_id, curPos++);
//       path_node.push_back(in_id);
//       // std::cerr << last_id << " " << in_id << "\n";
//       // abpoa_add_graph_edge(abg, last_id, in_id, 1, 1, add_read_id, 0, p_i, read_ids_n, p_n);
//       if (in_id < out_id) {
//         for (i = 0; i < out_id - in_id; ++i) {
//           // std::cerr << in_id + i << " " << in_id + i + 1 << "\n";
//           path_node.push_back(in_id + i + 1);
//           // DAG->add_adj(seqs.size(), in_id + i, in_id + i + 1, curPos++);
//         }
//       }
//       last_id = 0;
//       for (i = start_pos; i < end_pos; i++) {
//         DAG->add_adj(seqs.size(), last_id, path_node[i], curPos++);
//         last_id = path_node[i];
//       }
//       DAG->add_adj(seqs.size(), last_id, 1, curPos++);
//     }
//     // set abs
//     // abpoa_realloc_seq(abs);
//     seq_t tseq;
//     tseq.name = walk_name;
//     // abpoa_cpy_str(abs->name + abs->n_seq, walk_name, walk_name_len);
//     seqs.emplace_back(tseq);
//     // abs->is_rc[abs->n_seq] = is_rc; abs->n_seq++;
//   }
//   else err_fatal(__func__, "Error: no path in GFA path line (%s).", walk_name);
//   return 0;
// }
std::vector<seq_t> read_gfa(para_t *para, graph *DAG, const char *path, PathWriter * writer) {
  gzFile fp;
  kstring_t line = { 0,0,0 }, fa_seq = { 0,0,0 };
  kstream_t* ks;
  seg_seq_t* segs = seg_seq_init();
  int dret, is_fa = 0;
  // gfa_seg_t* fa_seg = 0;
  uint64_t lineno = 0;

  fp = path && strcmp(path, "-") ? gzopen(path, "r") : gzdopen(0, "r");
  std::vector<seq_t> seqs;
  if (fp == 0) return seqs;
  ks = ks_init(fp);
  // g = gfa_init();
  while (ks_getuntil(ks, KS_SEP_LINE, &line, &dret) >= 0) {
    // std::cerr << line.s << "\n";
    int ret = 0;
    ++lineno;
    if (line.l > 0 && line.s[0] == '>') { // FASTA header
      // is_fa = 1;
      // if (fa_seg) gfa_update_fa_seq(g, fa_seg, fa_seq.l, fa_seq.s);
      // fa_seg = gfa_parse_fa_hdr(g, s.s);
      // fa_seq.l = 0;
    }
    else if (is_fa) { // FASTA mode
      // if (s.l >= 3 && s.s[1] == '\t') { // likely a GFA line
      //   gfa_update_fa_seq(g, fa_seg, fa_seq.l, fa_seq.s); // finalize fa_seg
      //   fa_seg = 0;
      //   is_fa = 0;
      // }
      // else kputsn(s.s, s.l, &fa_seq); // likely a FASTA sequence line
    }
    if (is_fa) continue;
    if (line.l < 3 || line.s[1] != '\t') continue; // empty line
    if (line.s[0] == 'S') ret = gfa_parse_S(para, DAG, segs, line.s);
    // else if (s.s[0] == 'L') ret = gfa_parse_L(g, s.s);
    else if (line.s[0] == 'P') ret = gfa_parse_P(para, DAG, seqs, segs, line.s, writer);
    // else if (line.s[0] == 'W') ret = gfa_parse_W(para, DAG, seqs, segs, line.s);
    if (ret < 0)
      fprintf(stderr, "[E] invalid %c-line at line %ld (error code %d)\n", line.s[0], (long)lineno, ret);
  }
  // if (is_fa && fa_seg) gfa_update_fa_seq(g, fa_seg, fa_seq.l, fa_seq.s);
  // std::cerr << "topsort" << "\n";
  DAG->topsort(para, 0);
  free(fa_seq.s);
  free(line.s);
  // gfa_finalize(g);
  seg_seq_free(segs);
  ks_destroy(ks);
  gzclose(fp);
  return seqs;
}

// ============================================================================
// PathWriter Implementation
// ============================================================================

PathWriter::PathWriter(const char *filename) {
  unlink(filename);  // PathWrite clean
  fp_ = fopen(filename, "ab"); // Open in binary write mode
  if (!fp_) {
    fprintf(stderr, "[PathWriter] ERROR: Could not open file %s for writing.\n", filename);
  }
}

PathWriter::~PathWriter() {
  if (fp_) {
    fclose(fp_);
  }
}

void PathWriter::write_varint(uint32_t value) {
  unsigned char buffer[5];
  int i = 0;
  while (value >= 0x80) {
    buffer[i++] = (value & 0x7F) | 0x80;
    value >>= 7;
  }
  buffer[i++] = value;
  fwrite(buffer, 1, i, fp_);
}

bool PathWriter::write_path(int32_t seq_id, const std::vector<int> &node_ids) {
  if (!fp_) return false;

  // 1. Write Header: Sequence Order and Path Length
  uint32_t path_len = node_ids.size();
  if (fwrite(&seq_id, sizeof(int32_t), 1, fp_) != 1) return false;
  if (fwrite(&path_len, sizeof(uint32_t), 1, fp_) != 1) return false;

  // 2. Write Body: Node ID Deltas
  if (path_len == 0) return true;

  int last_node_id = 0;
  for (int current_node_id : node_ids) {
    int32_t delta = current_node_id - last_node_id;

    // 修复: 使用 ZigZag 编码将负数转为正数，避免 Varint 压缩失效占用 5 字节
    uint32_t zigzag_delta = (delta << 1) ^ (delta >> 31);

    write_varint(zigzag_delta);
    last_node_id = current_node_id;
  }
  return true;
}


// ============================================================================
// PathReader Implementation
// ============================================================================

PathReader::PathReader(const char *filename) {
  fp_ = fopen(filename, "rb"); // Open in binary read mode
  if (!fp_) {
    fprintf(stderr, "[PathReader] ERROR: Could not open file %s for reading.\n", filename);
  }
}

PathReader::~PathReader() {
  if (fp_) {
    fclose(fp_);
  }
}

bool PathReader::read_varint(uint32_t &value) {
  value = 0;
  int shift = 0;
  int c;

  // 修复: 使用 stdio 缓冲的 fgetc 替代单字节 fread，大幅提升读取性能
  while ((c = fgetc(fp_)) != EOF) {
    // 修复: 限制位移上限，防止恶意文件造成未定义行为和死循环
    if (shift >= 35) return false;

    value |= (uint32_t)(c & 0x7F) << shift;
    shift += 7;

    if ((c & 0x80) == 0) return true;
  }
  return false;
}

std::pair<int32_t, std::vector<int>> PathReader::read_next_path() {
  if (!fp_) {
    return { -1, {} };
  }

  // 预先嗅探是否到了文件末尾，避免干净的 EOF 引发错误日志
  int first_byte = fgetc(fp_);
  if (first_byte == EOF) {
    return { -1, {} };
  }
  ungetc(first_byte, fp_);

  // 1. Read Header
  int32_t seq_id; uint32_t path_len;
  if (fread(&seq_id, sizeof(int32_t), 1, fp_) != 1) {
    return { -1, {} };
  }
  if (fread(&path_len, sizeof(uint32_t), 1, fp_) != 1) {
    fprintf(stderr, "[PathReader] ERROR: Corrupted file - could not read path length.\n");
    return { seq_id, {} };
  }

  // 修复: 增加 path_len 的安全上限防护，防止解析错误长度直接造成 OOM 崩溃
  if (path_len > 100000000) {
    fprintf(stderr, "[PathReader] ERROR: Path length exceeds safety limits.\n");
    return { seq_id, {} };
  }

  // 2. Read Body
  std::vector<int> node_ids;
  node_ids.reserve(path_len);
  int last_node_id = 0;

  for (uint32_t i = 0; i < path_len; ++i) {
    uint32_t zigzag_delta;
    if (!read_varint(zigzag_delta)) {
      fprintf(stderr, "[PathReader] ERROR: Corrupted file - could not read varint delta.\n");
      return { seq_id, {} }; // Return empty vector on error
    }

    // 修复: ZigZag 解码，将无符号数还原为带符号的真实增量
    int32_t delta = (zigzag_delta >> 1) ^ -(int32_t)(zigzag_delta & 1);

    int current_node_id = last_node_id + delta;
    node_ids.push_back(current_node_id);
    last_node_id = current_node_id;
  }

  return { seq_id, node_ids };
}