clear 
close all
clc

mex -O -largeArrayDims -output gcd_mexed commonstuff.f90 gcd_mex.f90 gcd.f90 meshtopoaux.f90

