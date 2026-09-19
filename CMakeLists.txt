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
