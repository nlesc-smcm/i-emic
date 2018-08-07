rm -rf CMakeCache.txt CMakeFiles

cmake \
-DMRILU_DIR=/home/emulder/mrilu/ifort \
-DTrilinos_DIR=/home/emulder/trilinos/11.14 \
-DGTEST_DIR=/hpc/sw/gtest-1.7.0/ \
\
../
