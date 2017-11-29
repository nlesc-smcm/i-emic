rm -rf CMakeCache.txt CMakeFiles

cmake \
-DMRILU_DIR=/var/tmp/mrilu \
-DTrilinos_DIR=/var/tmp/trilinos/11.12/dynamic/lib/cmake/Trilinos \
\
../
