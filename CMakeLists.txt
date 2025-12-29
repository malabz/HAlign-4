cmake_minimum_required(VERSION 3.22)
project(halign4_viral)

set(CMAKE_CXX_STANDARD 17)

# 添加 WFA2 子模块（带有 CMake 支持）
add_subdirectory(PairwiseAlignment/WFA2-lib)

# 添加 include 路径（可选，wfa2cpp 会自动带上）
include_directories(
        PairwiseAlignment/WFA2-lib
)

# 添加你的可执行程序
add_executable(vhalign4
        main.cpp
        multi-thread/ThreadPoolPerThread.cpp
        multi-thread/ThreadPoolPerThread.hpp
        PairwiseAlignment/ViralPairwiseAlignment.cpp
        PairwiseAlignment/ViralPairwiseAlignment.hpp
        Utils/Arguments.cpp
        Utils/Arguments.hpp
        Utils/Fasta.cpp
        Utils/Fasta.hpp
        Utils/Graph.cpp
        Utils/Graph.hpp
        Utils/Insertion.cpp
        Utils/Insertion.hpp
        Utils/Utils.cpp
        Utils/Utils.hpp
        Utils/CommandLine.hpp
        Utils/Pseudo.hpp
)

# 正确链接 WFA2 的 C++ 静态库（推荐）
target_link_libraries(vhalign4 PRIVATE wfa2::wfa2cpp_static)
