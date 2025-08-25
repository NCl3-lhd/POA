cmake_minimum_required(VERSION 3.12)
set(CMAKE_CXX_STANDARD 17)  # 或 14、17、20

project(main)

# 设置可执行文件的输出目录为项目根目录下的 bin 文件夹
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# # SIMD
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")

# 1. 查找 zlib 包
find_package(ZLIB REQUIRED)

# 查找 src 目录下的所有 .cpp 文件
file(GLOB SOURCES "src/*.cpp" "include/*.c")

# 添加可执行文件
add_executable(POA ${SOURCES})

target_link_libraries(POA PRIVATE ZLIB::ZLIB)

# 设置头文件路径
target_include_directories(POA PRIVATE
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/src
)

# 设置默认构建类型为Release（自动包含-O3）
# 放在project()之后，任何目标定义之前
# if(NOT CMAKE_BUILD_TYPE)
#   set(CMAKE_BUILD_TYPE Release)
# endif()