CMAKE_MINIMUM_REQUIRED(VERSION 3.1)
project(opttp CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "-Ofast -march=native")

add_executable(opttp.out
  src/main.cpp
  )
