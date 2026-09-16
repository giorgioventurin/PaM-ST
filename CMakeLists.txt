cmake_minimum_required(VERSION 3.20)

project(pam_st CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Pin the optimisation level rather than relying on the compiler defaults CMake
# happens to supply. Debug and RelWithDebInfo keep their own levels, so a
# debugger or a sanitizer build still behaves.
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG" CACHE STRING "Release flags" FORCE)
endif()

find_package(Threads REQUIRED)

add_executable(pam_st
  src/main.cpp
  src/cli.cpp
  src/KDTree.cpp
  src/VPTreeJS.cpp
  src/utils.cpp
  src/io_utils.cpp
  src/null_model.cpp
  src/significance.cpp
)

target_include_directories(pam_st PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(pam_st PRIVATE Threads::Threads)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(pam_st PRIVATE -Wall -Wextra -Wpedantic)
endif()

include(CTest)
if(BUILD_TESTING)
  add_executable(block_null_unit
    tests/test_block_null.cpp
    src/null_model.cpp
  )
  target_include_directories(block_null_unit PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(block_null_unit PRIVATE -Wall -Wextra -Wpedantic)
  endif()
  add_test(NAME block_null_unit COMMAND block_null_unit)

  add_executable(significance_unit
    tests/test_significance.cpp
    src/significance.cpp
  )
  target_include_directories(significance_unit PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(significance_unit PRIVATE -Wall -Wextra -Wpedantic)
  endif()
  add_test(NAME significance_unit COMMAND significance_unit)
  add_test(
    NAME block_null_regression
    COMMAND ${CMAKE_COMMAND}
      -DPAM_ST=$<TARGET_FILE:pam_st>
      -DINPUT=${CMAKE_CURRENT_SOURCE_DIR}/tests/data/block_null_cells.csv
      -DFROZEN_INPUT=${CMAKE_CURRENT_SOURCE_DIR}/tests/data/frozen_cells.csv
      -DOUT_BASE=${CMAKE_CURRENT_BINARY_DIR}/block_null_regression
      -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_block_null.cmake
  )

  add_test(
    NAME frozen_population_smoke
    COMMAND pam_st
      --input ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/frozen_cells.csv
      --radius 1
      --rho 0
      --metric l2
      --permutations 5
      --max-motifs 2
      --freeze-cell-type P
  )
  set_tests_properties(
    frozen_population_smoke
    PROPERTIES PASS_REGULAR_EXPRESSION "candidate_centers=3"
  )

  add_test(
    NAME covering_neighborhood_smoke
    COMMAND pam_st
      --input ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/frozen_cells.csv
      --radius 1
      --rho 0
      --metric l2
      --permutations 5
      --max-motifs 2
      --neighborhood-mode covering
  )
  set_tests_properties(
    covering_neighborhood_smoke
    PROPERTIES PASS_REGULAR_EXPRESSION "coverage_fraction=1.000000"
  )

  add_test(
    NAME covering_freeze_smoke
    COMMAND pam_st
      --input ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/frozen_cells.csv
      --radius 1
      --rho 0
      --metric l2
      --permutations 5
      --max-motifs 2
      --neighborhood-mode covering
      --freeze-by-abundance
      --max-freeze-stages 1
  )
  set_tests_properties(
    covering_freeze_smoke
    PROPERTIES PASS_REGULAR_EXPRESSION "neighborhood_mode=covering"
  )

  add_test(
    NAME cumulative_freeze_sweep_smoke
    COMMAND pam_st
      --input ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/frozen_cells.csv
      --radius 1
      --rho 0
      --metric l2
      --permutations 5
      --max-motifs 2
      --freeze-by-abundance
      --stage-error-control holm
  )
  set_tests_properties(
    cumulative_freeze_sweep_smoke
    PROPERTIES PASS_REGULAR_EXPRESSION "freeze_stages=2"
  )

  add_test(
    NAME cumulative_freeze_regression
    COMMAND ${CMAKE_COMMAND}
      -DPAM_ST=$<TARGET_FILE:pam_st>
      -DINPUT=${CMAKE_CURRENT_SOURCE_DIR}/tests/data/frozen_cells.csv
      -DOUT_BASE=${CMAKE_CURRENT_BINARY_DIR}/freeze_sweep_regression
      -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_freeze_sweep.cmake
  )

  add_test(
    NAME normalized_anchor_regression
    COMMAND ${CMAKE_COMMAND}
      -DPAM_ST=$<TARGET_FILE:pam_st>
      -DINPUT=${CMAKE_CURRENT_SOURCE_DIR}/tests/data/proportional_histograms.csv
      -DOUT_BASE=${CMAKE_CURRENT_BINARY_DIR}/normalized_anchor_regression
      -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_normalized_anchors.cmake
  )
endif()
