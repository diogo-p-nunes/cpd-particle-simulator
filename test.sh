#!/bin/bash

make clean; make
echo "./simpar 1 3 1000000 20"
./simpar 1 3 1000000 20
#0.36 0.93
#0.50 0.50
echo "./simpar 1 10 2000000 10"
./simpar 1 10 2000000 10
#0.92 0.47
#0.50 0.50
echo "./simpar 1 30 20000000 10"
./simpar 1 30 20000000 10
#0.87 0.42
#0.50 0.50