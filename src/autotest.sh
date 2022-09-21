#!/bin/bash

clang++ -std=c++17 -O3 ./gen.cpp -o gen
Q=10
for i in `seq $Q`
do
  ./gen >$i.txt
done

clang++ -std=c++17 -O3 ./main.cpp -o main

for i in `seq $Q`
do
  ./main<$i.txt >out$i.txt &
done

wait

tail -n 1 -q out* > res
rm *.txt
