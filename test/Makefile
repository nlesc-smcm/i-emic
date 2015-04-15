include ../make.inc/Makefile.${PLAT}

internal_libs    = -l$(LIB_UTILS) -l$(LIB_OCEAN) -l$(LIB_NEWTON)
additional_flags = -Wall #--std=c++11
default: test

test: test.o 
	$(CXX) test.o $(CXX_FLAGS) -o test $(LIBRARY_DIRS) \
	$(internal_libs) $(LIBS)

%.o: %.C %.H
	$(CXX) -c $(additional_flags) $(INCLUDE_DIRS) $(CXX_FLAGS) $<
%.o: %.C
	$(CXX) -c $(additional_flags) $(INCLUDE_DIRS) $(CXX_FLAGS) $<
