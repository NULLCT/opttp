#!/bin/bash

clang++ -std=c++17 -O3 ./gen.cpp -o gen
./gen >1.txt
./gen >2.txt
./gen >3.txt
./gen >4.txt
./gen >5.txt
./gen >6.txt
./gen >7.txt
./gen >8.txt
./gen >9.txt
./gen >10.txt
./gen >11.txt
./gen >12.txt
./gen >13.txt
./gen >14.txt
./gen >15.txt
./gen >16.txt
./gen >17.txt
./gen >18.txt
./gen >19.txt
./gen >20.txt
./gen >21.txt
./gen >22.txt
./gen >23.txt
./gen >24.txt
./gen >25.txt

clang++ -std=c++17 -O3 ./main.cpp -o main

./main<1.txt >out1.txt &
./main<2.txt >out2.txt &
./main<3.txt >out3.txt &
./main<4.txt >out4.txt &
./main<5.txt >out5.txt &
./main<6.txt >out6.txt &
./main<7.txt >out7.txt &
./main<8.txt >out8.txt &
./main<9.txt >out9.txt &
./main<10.txt >out10.txt &
./main<11.txt >out11.txt &
./main<12.txt >out12.txt &
./main<13.txt >out13.txt &
./main<14.txt >out14.txt &
./main<15.txt >out15.txt &
./main<16.txt >out16.txt &
./main<17.txt >out17.txt &
./main<18.txt >out18.txt &
./main<19.txt >out19.txt &
./main<20.txt >out20.txt &
./main<21.txt >out21.txt &
./main<22.txt >out22.txt &
./main<23.txt >out23.txt &
./main<24.txt >out24.txt &
./main<25.txt >out25.txt &
