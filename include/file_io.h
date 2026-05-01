#ifndef FILE_IO_H
#define FILE_IO_H
#include <vector>
#include "sequence.h"
#include "parameter.h"
#include "graph.h"
#include <unistd.h>

class graph;
constexpr const char *tmp_path_file = "minipoa_paths.tmp";

void readFile(para_t* para, std::vector<seq_t>& seqs, const char* path);
void initPara(para_t* para);
std::vector<seq_t> read_gfa(para_t* para, graph* DAG, const char* path);
int gfa_parse_S(para_t* para, graph* DAG, seg_seq_t* segs, char* s);
int gfa_parse_W(para_t* para, graph* DAG, std::vector<seq_t>& seqs, seg_seq_t* segs, char* s);
int gfa_parse_P(para_t* para, graph* DAG, std::vector<seq_t>& seqs, seg_seq_t* segs, char* s);

struct PathWriter {

  PathWriter(const char *filename);
  ~PathWriter();

  // Writes the path for a single sequence to the file.
  // seq_id: The index of the sequence (0 to num_sequences - 1).
  // node_ids: A vector of node IDs that the sequence traverses.
  // Note: The caller is responsible for ensuring node_ids are in the correct traversal order.
  bool write_path(uint32_t seq_id, const std::vector<int> &node_ids);

  // Returns true if the writer was successfully opened.
  bool is_open() const { return fp_ != nullptr; }
  void close() { if (fp_) { fclose(fp_); fp_ = nullptr; } }

  FILE *fp_;

  // Helper function to encode an integer into varint format.
  void write_varint(uint32_t value);
};

struct PathReader {
  PathReader(const char *filename);
  ~PathReader();

  // Reads the path for the next sequence in the file.
  // Returns a pair containing the sequence order and its path of node IDs.
  // If end-of-file is reached or an error occurs, the vector will be empty.
  std::pair<uint32_t, std::vector<int>> read_next_path();

  // Returns true if the reader was successfully opened.
  bool is_open() const { return fp_ != nullptr; }

  FILE *fp_;

  // Helper function to decode a varint from the file stream.
  // Returns true on success, false on EOF or error.
  bool read_varint(uint32_t &value);
};

#endif