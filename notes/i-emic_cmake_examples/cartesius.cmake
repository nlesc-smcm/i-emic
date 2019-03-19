rm -rf CMakeCache.txt CMakeFiles

cmake \
-D CMAKE_C_COMPILER:STRING="mpiicc" \
-D CMAKE_CXX_COMPILER:STRING="mpiicpc" \
-D CMAKE_Fortran_COMPILER:STRING="mpiifort" \
-D Trilinos_DIR=/${HOME}/trilinos/12.12  \
../
