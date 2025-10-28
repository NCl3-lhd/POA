#include <zlib.h>
#include <stdexcept> 
#include <stdint.h>
#include "handle_input.h"
#include "kseq.h"
#include "sequence.h" 
#include "utils.h"

KSEQ_INIT(gzFile, gzread)

// handle input file
void readFile(std::vector<seq_t>& seqs, const char* path) {
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
int gfa_parse_W(para_t* para, graph* DAG, std::vector<seq_t>& seqs, seg_seq_t* segs, char* s) {
  if (s[1] != '\t' || s[2] == '\0') return -1;
  char* end_s, * start_s, * path = 0;
  int i, is_ok = 0, is_rc = -1;
  char* walk_name = 0; int walk_name_len = 0;
  kstring_t* seg_seq, * seg_name;
  // std::cerr << s << "\n";
  for (i = 0, end_s = start_s = s + 2;; ++end_s) {
    if (*end_s == '\0' || *end_s == '\t') {
      int c = *end_s;
      *end_s = '\0';
      if (i == 0) {
        // walk id
      }
      else if (i == 1) {

      }
      else if (i == 2) {
        walk_name = start_s;
        walk_name_len = end_s - start_s;
      }
      else if (i == 3) {

      }
      else if (i == 4) {

      }
      else if (i == 5) {
        path = start_s;
        is_ok = 1;
        break;
      }
      if (c == 0) break;
      ++i, start_s = end_s + 1;
    }
  }
  // std::cerr << walk_name << " " << path << "\n";
  // std::cerr << "ok" << "\n";
  if (is_ok) {
    char* end_s, * start_s, * _seg_name; khint_t pos, seg_pos; int absent;
    int id, in_id = -1, out_id = -1, last_id = 0, next_id = 1;
    int curPos = 0;
    for (end_s = start_s = path + 1; ; ++end_s) {
      if (*end_s == '>') {
        if (is_rc == 1) err_fatal(__func__, "Error: walk has both \'>\' and \'<\' seg. (%s)", walk_name);
        is_rc = 0; *end_s = '\0'; _seg_name = start_s;
        // std::cerr << _seg_name << " ";
        // if (_seg_name == "s211") {
        //   exit(1);
        // }
        seg_pos = kh_get(str, segs->h, _seg_name);
        if (seg_pos == kh_end(segs->h)) err_fatal(__func__, "Error: seg (%s) not exist.", start_s);
        seg_name = segs->name + kh_val(segs->h, seg_pos);
        seg_seq = segs->seq + kh_val(segs->h, seg_pos);

        // check if seg already exist
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
        //   in_id = kh_val(seg_name2in_id, pos);
        //   pos = kh_put(str, seg_name2out_id, seg_name->s, &absent);
        //   out_id = kh_val(seg_name2out_id, pos);
        // }
        // add edge
        // std::cerr << in_id << " " << out_id << "\n";
        DAG->add_adj(seqs.size(), last_id, in_id, curPos++);
        // std::cerr << last_id << " " << in_id << "\n";
        // abpoa_add_graph_edge(abg, last_id, in_id, 1, 1, add_read_id, 0, p_i, read_ids_n, p_n);
        if (in_id < out_id) {
          for (i = 0; i < out_id - in_id; ++i) {
            // std::cerr << in_id + i << " " << in_id + i + 1 << "\n";
            DAG->add_adj(seqs.size(), in_id + i, in_id + i + 1, curPos++);
          }
        }
        else if (in_id > out_id) err_fatal(__func__, "Error: in_id (%d) > out_id (%d).", in_id, out_id);

        last_id = out_id;
        start_s = end_s + 1;
      }
      else if (*end_s == '<') { // is not correct, not support <
        if (is_rc == 0) err_fatal(__func__, "Error: walk has both \'>\' and \'<\' seg. (%s)", walk_name);
        is_rc = 1; *end_s = '\0'; _seg_name = start_s;
        seg_pos = kh_get(str, segs->h, _seg_name);
        if (seg_pos == kh_end(segs->h)) err_fatal(__func__, "Error: seg (%s) not exist.", start_s);
        seg_name = segs->name + kh_val(segs->h, seg_pos);
        seg_seq = segs->seq + kh_val(segs->h, seg_pos);

        // check if seg exist
        in_id = segs->start_id.a[seg_pos];
        out_id = segs->end_id.a[seg_pos];
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
        DAG->add_adj(seqs.size(), out_id, next_id, curPos++);
        // abpoa_add_graph_edge(abg, out_id, next_id, 1, 1, add_read_id, 0, p_i, read_ids_n, p_n);
        if (in_id < out_id) {
          for (i = 0; i < out_id - in_id; ++i)
            DAG->add_adj(seqs.size(), in_id + i, in_id + i + 1, curPos++);
          // abpoa_add_graph_edge(abg, in_id + i, in_id + i + 1, 1, 1, add_read_id, 0, p_i, read_ids_n, p_n);
        }
        else if (in_id > out_id) err_fatal(__func__, "Error: in_id (%d) > out_id (%d).", in_id, out_id);

        next_id = in_id;
        start_s = end_s + 1;
      }
      else if (*end_s == '\0' || *end_s == '\t') break;
    }

    if (is_rc) {
      // abpoa_add_graph_edge(abg, ABPOA_SRC_NODE_ID, next_id, 1, 1, add_read_id, 0, p_i, read_ids_n, p_n);
      DAG->add_adj(seqs.size(), 0, next_id, curPos++);
    }
    else {
      //abpoa_add_graph_edge(abg, last_id, ABPOA_SINK_NODE_ID, 1, 1, add_read_id, 0, p_i, read_ids_n, p_n);
      // std::cerr << start_s << " ";
      seg_pos = kh_get(str, segs->h, start_s);
      if (seg_pos == kh_end(segs->h)) err_fatal(__func__, "Error: seg (%s) not exist.", start_s);
      seg_seq = segs->seq + kh_val(segs->h, seg_pos);

      // check if seg already exist
      in_id = segs->start_id.a[kh_val(segs->h, seg_pos)];
      out_id = segs->end_id.a[kh_val(segs->h, seg_pos)];
      DAG->add_adj(seqs.size(), last_id, in_id, curPos++);
      // std::cerr << last_id << " " << in_id << "\n";
      // abpoa_add_graph_edge(abg, last_id, in_id, 1, 1, add_read_id, 0, p_i, read_ids_n, p_n);
      if (in_id < out_id) {
        for (i = 0; i < out_id - in_id; ++i) {
          // std::cerr << in_id + i << " " << in_id + i + 1 << "\n";
          DAG->add_adj(seqs.size(), in_id + i, in_id + i + 1, curPos++);
        }
      }
      last_id = out_id;
      DAG->add_adj(seqs.size(), last_id, 1, curPos++);
    }
    // set abs
    // abpoa_realloc_seq(abs);
    seq_t tseq;
    tseq.name = walk_name;
    // abpoa_cpy_str(abs->name + abs->n_seq, walk_name, walk_name_len);
    seqs.emplace_back(tseq);
    // abs->is_rc[abs->n_seq] = is_rc; abs->n_seq++;
  }
  else err_fatal(__func__, "Error: no path in GFA path line (%s).", walk_name);
  return 0;
}
std::vector<seq_t> read_gfa(para_t* para, graph* DAG, const char* path) {
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
    else if (line.s[0] == 'W') ret = gfa_parse_W(para, DAG, seqs, segs, line.s);
    if (ret < 0)
      fprintf(stderr, "[E] invalid %c-line at line %ld (error code %d)\n", line.s[0], (long)lineno, ret);
  }
  // if (is_fa && fa_seg) gfa_update_fa_seq(g, fa_seg, fa_seq.l, fa_seq.s);
  // std::cerr << "topsort" << "\n";
  DAG->topsort(0, para->f);
  free(fa_seq.s);
  free(line.s);
  // gfa_finalize(g);
  seg_seq_free(segs);
  ks_destroy(ks);
  gzclose(fp);
  return seqs;
}