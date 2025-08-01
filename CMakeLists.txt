cmake_minimum_required(VERSION 3.12)

project(main)

# 查找 src 目录下的所有 .cpp 文件
file(GLOB SOURCES "src/*.cpp")
find_package(ZLIB REQUIRED)

# 添加可执行文件
add_executable(POA ${SOURCES})

# 设置头文件路径
target_include_directories(POA PRIVATE
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/src
)

# 设置默认构建类型为Release（自动包含-O3）
# 放在project()之后，任何目标定义之前
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()