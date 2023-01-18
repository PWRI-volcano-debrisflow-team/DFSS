del DR_ver_1.0.f90
del DR.exe
copy ..\..\..\DFSS_source_code\DR_ver_1.0.f90 .
gfortran -o DR.exe -fopenmp DR_ver_1.0.f90
set NUM_OMP_THREADS=8
DR.exe
