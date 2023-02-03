CMAKE_MINIMUM_REQUIRED(VERSION 3.1)
project(opttp CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "-Ofast -march=native")

add_executable(opttp.out
  src/main.cpp
  src/directed_graph.cpp
  src/simulated_annealing.cpp
  src/unionfind.cpp
  )

enable_testing()
add_executable(test_submodules
  src/test.cpp
  src/directed_graph.cpp
  src/simulated_annealing.cpp
  src/unionfind.cpp
)

add_test(NAME test_sub COMMAND test_submodules)