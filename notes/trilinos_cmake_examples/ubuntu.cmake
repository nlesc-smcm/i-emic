rm -rf CMakeCache.txt CMakeFiles

cmake \
   -D CMAKE_INSTALL_PREFIX:PATH=/home/erik/local \
   -D TPL_ENABLE_MPI:BOOL=ON \
   -D TPL_ENABLE_METIS:BOOL=ON \
   -D TPL_ENABLE_ParMETIS:BOOL=ON \
   -D TPL_ENABLE_HDF5:BOOL=ON \
\
   -D Trilinos_VERBOSE_CONFIGURE:BOOL=OFF \
   -D CMAKE_VERBOSE_MAKEFILE:BOOL=ON \
\
   -D CMAKE_BUILD_TYPE:STRING=RELEASE \
   -D BUILD_SHARED_LIBS=ON \
   -D CMAKE_CXX_FLAGS:STRING="-g -O3" \
   -D Trilinos_ENABLE_DEBUG=OFF \
\
   -D Trilinos_ENABLE_Export_Makefiles:BOOL=ON \
   -D Trilinos_ENABLE_EXAMPLES:BOOL=OFF \
   -D Trilinos_ENABLE_TESTS:BOOL=OFF \
   -D Belos_ENABLE_Experimental:BOOL=OFF \
\
   -D Trilinos_ENABLE_Teuchos:BOOL=ON \
   -D Trilinos_ENABLE_Epetra:BOOL=ON \
   -D Trilinos_ENABLE_EpetraExt:BOOL=ON \
   -D Trilinos_ENABLE_Ifpack:BOOL=ON \
   -D Trilinos_ENABLE_Amesos:BOOL=ON \
   -D Trilinos_ENABLE_Anasazi:BOOL=ON \
   -D Trilinos_ENABLE_Belos:BOOL=ON \
   -D Trilinos_ENABLE_ML:BOOL=ON \
   -D Trilinos_ENABLE_OpenMP:BOOL=OFF \
\
../
