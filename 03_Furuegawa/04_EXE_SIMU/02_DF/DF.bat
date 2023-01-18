del DF_ver_1.0.f90
del DF.exe
copy ..\..\..\DFSS_source_code\DF_ver_1.0.f90 .
gfortran -o DF.exe -fopenmp DF_ver_1.0.f90
set NUM_OMP_THREADS=8
DF.exe
