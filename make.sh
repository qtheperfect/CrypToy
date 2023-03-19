#! /bin/bash
echo "Requiring: libgmp-dev, libgmp10, libgmpxx4ldbl"
echo ; echo ; 
cxxfile="main.cpp"
g++ -O3 -o main.exe -std=c++2a -fconcepts "$cxxfile" -lgmp -lgmpxx && ./main.exe
echo; echo;
