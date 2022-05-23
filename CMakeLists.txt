CMAKE_MINIMUM_REQUIRED(VERSION 3.14)
project(opttp CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -fsanitize=undefined,address -g -O3")

add_executable(opttp.out
  src/main.cpp
  )
