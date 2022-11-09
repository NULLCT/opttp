CMAKE_MINIMUM_REQUIRED(VERSION 3.24)
project(opttp CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "-Wall -Ofast")

add_executable(opttp.out
  src/main.cpp
  )
