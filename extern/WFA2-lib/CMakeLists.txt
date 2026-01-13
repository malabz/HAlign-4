# For Debian currently with
#
#   cd build
#   cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo -DOPENMP=1 ..
#   cmake --build .
#   ctest .
#   make install
# See below option statements and the README for build information

cmake_minimum_required(VERSION 3.16)
project(wfa2lib)

set(CMAKE_CXX_STANDARD 17)

include(FeatureSummary)
include(GNUInstallDirs)

# ---------------------------------------------------------------------------
# 作为子项目被 add_subdirectory() 引入时，主工程通常只需要核心库。
# 但 WFA2-lib 上游 CMake 默认总是添加 benchmark 可执行文件，并注册测试。
#
# 这里补充两个构建开关：
# - WFA2LIB_BUILD_BENCHMARK：是否构建 align_benchmark（默认 ON，保持上游行为不变）
# - WFA2LIB_BUILD_TESTS：是否 enable_testing 并 add_test（默认 ON，保持上游行为不变）
#
# 主工程可以在 add_subdirectory 之前 set(... CACHE BOOL ... FORCE) 来关闭它们，
# 从而“只编译核心代码”。
# ---------------------------------------------------------------------------
option(WFA2LIB_BUILD_BENCHMARK "Build WFA2-lib benchmark executable" ON)
option(WFA2LIB_BUILD_TESTS "Build/register WFA2-lib tests" ON)

# ---------------------------------------------------------------------------
# PkgConfig 并不是 WFA2-lib 核心库（wfa2/wfa2cpp）的必需依赖。
#
# 你的主工程希望完全不使用 pkg-config（例如设置了
# CMAKE_DISABLE_FIND_PACKAGE_PkgConfig=ON），因此这里不能用 REQUIRED，
# 否则会导致 CMake 配置阶段直接失败。
#
# 改为 QUIET：
# - 有 pkg-config：照常检测（保持与上游尽量一致）
# - 没有/被禁用：继续编译核心库，仅打印提示
# ---------------------------------------------------------------------------
find_package(PkgConfig QUIET)
if(NOT PkgConfig_FOUND)
  message(STATUS "PkgConfig not found (or disabled). Proceeding without it.")
endif()

feature_summary(
  FATAL_ON_MISSING_REQUIRED_PACKAGES
  WHAT REQUIRED_PACKAGES_NOT_FOUND)

# ---- Options

option(OPENMP "Enable OpenMP" OFF) # enables WFA_PARALLEL
option(PROFILING "Enable profiling" OFF)
option(ASAN "Use address sanitiser" OFF)
option(EXTRA_FLAGS "Add optimization flags for C/C++ compiler" OFF)

# include(CheckIPOSupported) # adds lto
# check_ipo_supported(RESULT ipo_supported OUTPUT output)

# ---- Dependencies

if(OPENMP)
  include(FindOpenMP)
  set(OPTIMIZE_FLAGS "-DWFA_PARALLEL")
endif(OPENMP)

if(EXTRA_FLAGS)
  set(OPTIMIZE_FLAGS "${OPTIMIZE_FLAGS} ${EXTRA_FLAGS}")
endif(EXTRA_FLAGS)

find_package(Threads)
set_package_properties(Threads PROPERTIES TYPE REQUIRED)

# ---- Build switches
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
# set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ${ipo_supported})

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING
          "Choose the type of build, options are: Release|Debug|RelWithDebInfo (for distros)." FORCE)
endif()

if (${CMAKE_BUILD_TYPE} MATCHES Release)
  set(OPTIMIZE_FLAGS "${OPTIMIZE_FLAGS} -march=native -D_FILE_OFFSET_BITS=64")
endif()

if ((${CMAKE_BUILD_TYPE} MATCHES Release) OR (${CMAKE_BUILD_TYPE} MATCHES RelWithDebInfo))
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OPTIMIZE_FLAGS}")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OPTIMIZE_FLAGS}")
endif ()

if (${CMAKE_BUILD_TYPE} MATCHES "Debug")
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OPTIMIZE_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OPTIMIZE_FLAGS}")
  add_definitions(-Wfatal-errors)
endif ()

if (ASAN)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=address  -fno-omit-frame-pointer -fno-common")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address  -fno-omit-frame-pointer -fno-common")
endif(ASAN)

if(PROFILING)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
endif(PROFILING)

if(GPROF)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
endif(GPROF)

# ---- Include files

file(GLOB INCLUDES
  wavefront/*.h*
  )
file(GLOB UTILS_INCLUDES
  utils/*.h*
  )
file(GLOB ALIGNMENT_INCLUDES
  alignment/*.h*
  )
file(GLOB SYSTEM_INCLUDES
  system/*.h*
  )

set(wfa2lib_SOURCE
  wavefront/wavefront_align.c
  wavefront/wavefront_aligner.c
  wavefront/wavefront_attributes.c
  wavefront/wavefront_backtrace_buffer.c
  wavefront/wavefront_backtrace.c
  wavefront/wavefront_backtrace_offload.c
  wavefront/wavefront_bialign.c
  wavefront/wavefront_bialigner.c
  wavefront/wavefront.c
  wavefront/wavefront_components.c
  wavefront/wavefront_compute_affine2p.c
  wavefront/wavefront_compute_affine.c
  wavefront/wavefront_compute.c
  wavefront/wavefront_compute_edit.c
  wavefront/wavefront_compute_linear.c
  wavefront/wavefront_debug.c
  wavefront/wavefront_display.c
  wavefront/wavefront_extend.c
  wavefront/wavefront_heuristic.c
  wavefront/wavefront_pcigar.c
  wavefront/wavefront_penalties.c
  wavefront/wavefront_plot.c
  wavefront/wavefront_sequences.c
  wavefront/wavefront_slab.c
  wavefront/wavefront_unialign.c
  wavefront/wavefront_termination.c
  wavefront/wavefront_extend_kernels_avx.c
  wavefront/wavefront_extend_kernels.c
  system/mm_stack.c
  system/mm_allocator.c
  system/profiler_counter.c
  system/profiler_timer.c
  utils/bitmap.c
  utils/dna_text.c
  utils/sequence_buffer.c
  utils/vector.c
  utils/commons.c
  utils/heatmap.c
  alignment/affine2p_penalties.c
  alignment/affine_penalties.c
  alignment/cigar.c
  alignment/cigar_utils.c
  alignment/score_matrix.c
)

add_library(wfa2_static
    ${wfa2lib_SOURCE}
    )
add_library(wfa2 SHARED ${wfa2lib_SOURCE})
set_target_properties(wfa2_static PROPERTIES OUTPUT_NAME wfa2)
set_target_properties(wfa2 PROPERTIES SOVERSION 0)
target_include_directories(wfa2 PUBLIC . wavefront utils)
target_include_directories(wfa2_static PUBLIC . wavefront utils)
add_library(wfa2::wfa2 ALIAS wfa2)
add_library(wfa2::wfa2_static ALIAS wfa2_static)

if(OPENMP)
  target_link_libraries(wfa2_static PRIVATE OpenMP::OpenMP_C)
  target_link_libraries(wfa2 PRIVATE OpenMP::OpenMP_C)
endif(OPENMP)

# ---- C++ binding library

set(wfa2cpp_SOURCE
  bindings/cpp/WFAligner.cpp
)
file(GLOB CPP_INCLUDES
  bindings/cpp/*.h*
  )
add_library(wfa2cpp_static STATIC ${wfa2cpp_SOURCE})
add_library(wfa2cpp SHARED ${wfa2cpp_SOURCE})
set_target_properties(wfa2cpp PROPERTIES SOVERSION 0)
set_target_properties(wfa2cpp_static PROPERTIES OUTPUT_NAME wfa2cpp)
target_link_libraries(wfa2cpp PUBLIC wfa2)
target_link_libraries(wfa2cpp_static PUBLIC wfa2_static)
add_library(wfa2::wfa2cpp ALIAS wfa2cpp)
add_library(wfa2::wfa2cpp_static ALIAS wfa2cpp_static)

if(OPENMP)
  target_link_libraries(wfa2cpp_static PRIVATE OpenMP::OpenMP_CXX)
  target_link_libraries(wfa2cpp PRIVATE OpenMP::OpenMP_CXX)
endif(OPENMP)

# ---- Get version

file (STRINGS "VERSION.txt" BUILD_NUMBER)
add_definitions(-DWFA2LIB_VERSION="${BUILD_NUMBER}")
add_definitions(-DVERSION="${BUILD_NUMBER}")

set(wfa2lib_LIBS
)

# add_dependencies(wfa2lib ${wfa2lib_DEPS})

# ---- Build all
# add_subdirectory(tools/align_benchmark)

if(WFA2LIB_BUILD_BENCHMARK)
  add_executable(align_benchmark
    tools/align_benchmark/align_benchmark_params.c
    tools/align_benchmark/align_benchmark.c
    tools/align_benchmark/benchmark/benchmark_check.c
    tools/align_benchmark/benchmark/benchmark_edit.c
    tools/align_benchmark/benchmark/benchmark_gap_affine2p.c
    tools/align_benchmark/benchmark/benchmark_gap_affine.c
    tools/align_benchmark/benchmark/benchmark_gap_linear.c
    tools/align_benchmark/benchmark/benchmark_indel.c
    tools/align_benchmark/benchmark/benchmark_utils.c
    system/mm_allocator.c
    system/mm_stack.c
    system/profiler_counter.c
    system/profiler_timer.c
    utils/bitmap.c
    utils/commons.c
    utils/dna_text.c
    utils/heatmap.c
    utils/sequence_buffer.c
    utils/vector.c
    alignment/affine2p_penalties.c
    alignment/affine_penalties.c
    alignment/cigar.c
    alignment/score_matrix.c
    tools/align_benchmark/gap_affine/affine_matrix.c
    tools/align_benchmark/gap_affine/swg.c
    tools/align_benchmark/gap_linear/nw.c
    tools/align_benchmark/gap_affine2p/affine2p_dp.c
    tools/align_benchmark/gap_affine2p/affine2p_matrix.c
    tools/align_benchmark/edit/edit_bpm.c
    tools/align_benchmark/edit/edit_dp.c
    tools/align_benchmark/indel/indel_dp.c
  )
  target_include_directories(align_benchmark PUBLIC . wavefront system utils tools/align_benchmark)
  add_dependencies(align_benchmark wfa2)
  target_link_libraries(align_benchmark m)
  target_link_libraries(align_benchmark wfa2)
  if(OPENMP)
    target_link_libraries(align_benchmark OpenMP::OpenMP_C)
  endif(OPENMP)
endif()

# ---- Test

if(WFA2LIB_BUILD_TESTS)
  enable_testing()

  function(add_wfa_test)
    add_test(
      NAME wfa2lib
      COMMAND ./tests/wfa.utest.sh ${CMAKE_CURRENT_BINARY_DIR}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      )
  endfunction()

  add_wfa_test()
endif()

# ---- Install

# Do not install anything when used with FetchContent
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  install(TARGETS wfa2_static ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} )
  install(TARGETS wfa2 ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} )

  install(FILES ${INCLUDES} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/wfa2lib/wavefront)

  install(FILES ${UTILS_INCLUDES} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/wfa2lib/utils)

  install(FILES ${ALIGNMENT_INCLUDES} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/wfa2lib/alignment)

  install(FILES ${SYSTEM_INCLUDES} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/wfa2lib/system)

  install(FILES ${CPP_INCLUDES} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/wfa2lib/cpp)

  install(TARGETS wfa2cpp_static ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} )
  install(TARGETS wfa2cpp ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBdir} )
endif()
