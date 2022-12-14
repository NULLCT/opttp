CMAKE_MINIMUM_REQUIRED(VERSION 3.24)
project(opttp CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS_DEBUG "-Og -g -D DEBUG -fsanitize=address,undefined,float-divide-by-zero,bounds,float-cast-overflow,alignment,object-size,vptr,bounds-strict,signed-integer-overflow,builtin,pointer-overflow,pointer-compare,pointer-subtract")
set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -march=native")

add_executable(opttp.out
  src/main.cpp
  )
add_executable(gen.out
  src/gen.cpp
  )
