CMAKE_MINIMUM_REQUIRED(VERSION 3.14)
project(opttp CXX)

set(CMAKE_CXX_STANDARD 17)
if(UNIX)
  set(CMAKE_CXX_FLAGS "-Wall -Wextra -fsanitize=undefined,address,memory -g")
endif(UNIX)
if(APPLE)
  set(CMAKE_CXX_FLAGS "-Wall -Wextra -fsanitize=undefined,address -g")
endif(APPLE)

add_executable(opttp.out
  src/main.cpp
  )
