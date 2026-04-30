#ifndef SIMDPRIV_H
#define SIMDPRIV_H
#include <immintrin.h>
#include <cstddef>

// ============================================================
//  SIMD backend selection (controlled by CMake)
// ============================================================

// ============================================================
//  SIMD register & width
// ============================================================

#if defined(ENABLE_AVX512)
using simd_reg = __m512i;
constexpr size_t SIMD_BYTES = 64;
#elif defined(ENABLE_AVX2)
using simd_reg = __m256i;
constexpr size_t SIMD_BYTES = 32;
#else
using simd_reg = __m128i;
constexpr size_t SIMD_BYTES = 16;
#endif

constexpr int simd_width = SIMD_BYTES / sizeof(int);

// ============================================================
//  Load / Store
// ============================================================

inline simd_reg simd_load(const int* p) {
#if defined(ENABLE_AVX512)
  return _mm512_load_si512(p);
#elif defined(ENABLE_AVX2)
  return _mm256_load_si256(reinterpret_cast<const __m256i*>(p));
#else
  return _mm_load_si128(reinterpret_cast<const __m128i*>(p));
#endif
}

inline simd_reg simd_loadu(const int* p) {
#if defined(ENABLE_AVX512)
  return _mm512_loadu_si512(p);
#elif defined(ENABLE_AVX2)
  return _mm256_loadu_si256(reinterpret_cast<const __m256i*>(p));
#else
  return _mm_loadu_si128(reinterpret_cast<const __m128i*>(p));
#endif
}

inline void simd_store(int* p, simd_reg v) {
#if defined(ENABLE_AVX512)
  _mm512_store_si512(p, v);
#elif defined(ENABLE_AVX2)
  _mm256_store_si256(reinterpret_cast<__m256i*>(p), v);
#else
  _mm_store_si128(reinterpret_cast<__m128i*>(p), v);
#endif
}

// ============================================================
//  Set
// ============================================================

inline simd_reg simd_set1(int x) {
#if defined(ENABLE_AVX512)
  return _mm512_set1_epi32(x);
#elif defined(ENABLE_AVX2)
  return _mm256_set1_epi32(x);
#else
  return _mm_set1_epi32(x);
#endif
}

// ============================================================
//  Arithmetic
// ============================================================

inline simd_reg simd_add(simd_reg a, simd_reg b) {
#if defined(ENABLE_AVX512)
  return _mm512_add_epi32(a, b);
#elif defined(ENABLE_AVX2)
  return _mm256_add_epi32(a, b);
#else
  return _mm_add_epi32(a, b);
#endif
}

inline simd_reg simd_max(simd_reg a, simd_reg b) {
#if defined(ENABLE_AVX512)
  return _mm512_max_epi32(a, b);
#elif defined(ENABLE_AVX2)
  return _mm256_max_epi32(a, b);
#else
  return _mm_max_epi32(a, b);
#endif
}

// ============================================================
//  DP helper
//  lane 0 = prev
//  lane i = p[i-1]
// ============================================================

inline simd_reg simd_set_prev_and_load(int prev, const int* p) {
#if defined(ENABLE_AVX512)
  return _mm512_setr_epi32(
      prev,
      p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7],
      p[8], p[9], p[10], p[11], p[12], p[13], p[14]
  );
#elif defined(ENABLE_AVX2)
  return _mm256_setr_epi32(
      prev,
      p[0], p[1], p[2], p[3], p[4], p[5], p[6]
  );
#else
  return _mm_setr_epi32(
      prev,
      p[0], p[1], p[2]
  );
#endif
}

// ============================================================
//  ArgMax (Find index of maximum value)
// ============================================================

// 辅助函数：初始化 0 到 simd_width-1 的递增下标序列
inline simd_reg simd_set_sequence() {
#if defined(ENABLE_AVX512)
  return _mm512_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15);
#elif defined(ENABLE_AVX2)
  return _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7);
#else
  return _mm_setr_epi32(0, 1, 2, 3);
#endif
}

// 核心函数：向量化寻找最大值的下标
inline int simd_argmax(const int* M, int block_num) {
  if (block_num <= 0) return 0;

  // 初始化：加载第一个块的值，并设定初始对应的下标
  simd_reg max_vals = simd_load(M);
  simd_reg max_idxs = simd_set_sequence();
  
  simd_reg curr_idxs = max_idxs;
  simd_reg idx_inc = simd_set1(simd_width);

  // 遍历后续所有的 SIMD 块
  for (int i = 1; i < block_num; ++i) {
    curr_idxs = simd_add(curr_idxs, idx_inc); // 更新当前块的真实基础下标
    simd_reg curr_vals = simd_load(M + i * simd_width);

#if defined(ENABLE_AVX512)
    // AVX-512: 使用 Mask 寄存器进行高效率的 Blend
    __mmask16 mask = _mm512_cmpgt_epi32_mask(curr_vals, max_vals);
    max_vals = _mm512_mask_blend_epi32(mask, max_vals, curr_vals);
    max_idxs = _mm512_mask_blend_epi32(mask, max_idxs, curr_idxs);
#elif defined(ENABLE_AVX2)
    // AVX2: 使用 _mm256_blendv_epi8，利用 _mm256_cmpgt_epi32 生成的 32-bit 掩码进行条件赋值
    __m256i mask = _mm256_cmpgt_epi32(curr_vals, max_vals);
    max_vals = _mm256_blendv_epi8(max_vals, curr_vals, mask);
    max_idxs = _mm256_blendv_epi8(max_idxs, curr_idxs, mask);
#else
    // SSE2: 没有原生按位 Blend 指令，使用标准的位运算 (AND/ANDNOT/OR) 组合实现条件替换
    __m128i mask = _mm_cmpgt_epi32(curr_vals, max_vals);
    max_vals = _mm_or_si128(_mm_and_si128(mask, curr_vals),
                            _mm_andnot_si128(mask, max_vals));
    max_idxs = _mm_or_si128(_mm_and_si128(mask, curr_idxs),
                            _mm_andnot_si128(mask, max_idxs));
#endif
  }

  // 收尾阶段：水平归约 (Horizontal Reduction)
  // 将包含局部最大值及其下标的 SIMD 寄存器存入内存，进行小规模的标量归约
  alignas(SIMD_BYTES) int final_vals[simd_width];
  alignas(SIMD_BYTES) int final_idxs[simd_width];
  simd_store(final_vals, max_vals);
  simd_store(final_idxs, max_idxs);

  int best_idx = final_idxs[0];
  int best_val = final_vals[0];

  // 在 simd_width (通常为 4, 8 或 16) 个元素中找出真正的全局最大值
  for (int i = 1; i < simd_width; ++i) {
    // 处理逻辑：发现更大的值，或者值相等但下标更靠前的情况（完美还原标量版本的逻辑行为）
    if (final_vals[i] > best_val || 
       (final_vals[i] == best_val && final_idxs[i] < best_idx)) {
      best_val = final_vals[i];
      best_idx = final_idxs[i];
    }
  }

  return best_idx;
}

#endif
