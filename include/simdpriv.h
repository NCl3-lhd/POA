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
#endif