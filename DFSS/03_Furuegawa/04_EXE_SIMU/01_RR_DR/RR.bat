del RR_ver_1.0.f90
del RR.exe
copy ..\..\..\DFSS_source_code\RR_ver_1.0.f90 .
gfortran -o RR.exe -fopenmp RR_ver_1.0.f90
set NUM_OMP_THREADS=8
RR.exe
