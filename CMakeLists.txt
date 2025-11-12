cmake_minimum_required(VERSION 3.12)
set(CMAKE_CXX_STANDARD 17)  # 或 14、17、20

project(main)

# 设置可执行文件的输出目录为项目根目录下的 bin 文件夹
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# # SIMD
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")

# 1. 查找 zlib 包
find_package(ZLIB REQUIRED)

find_package(OpenMP)
add_compile_options(-Wunknown-pragmas)


# 查找 src 目录下的所有 .cpp 文件
file(GLOB SOURCES "src/*.cpp" "src/*.c")

# 添加可执行文件
add_executable(minipoa ${SOURCES})

target_link_libraries(minipoa PRIVATE ZLIB::ZLIB)
target_link_libraries(minipoa PRIVATE OpenMP::OpenMP_CXX)
# 设置头文件路径
target_include_directories(minipoa PRIVATE
    ${PROJECT_SOURCE_DIR}/include
)

# 设置默认构建类型为Release（自动包含-O3）
# 放在project()之后，任何目标定义之前
# if(NOT CMAKE_BUILD_TYPE)
#   set(CMAKE_BUILD_TYPE Release)
# endif()