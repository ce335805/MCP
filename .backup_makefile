CC = icpc
OPT_REPORT = -qopt-report=4 -qopt-report-phase ipo
CFLAGS = -std=c++11 -O3 -Wall -Wextra
MKL_SEQ = -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -liomp5 -lpthread
OBJ_FILES = utils.o costs.o MY_MCP.o
SRC_FILES = utils.cpp costs.cpp MY_MCP.cpp

MY_MCP: $(OBJ_FILES)
	$(CC) $(CFLAGS) $(MKL_SEQ) $(OBJ_FILES) -o $@

MY_MCP.o: MY_MCP.cpp costs.cpp costs.h utils.cpp utils.h
	$(CC) $(CFLAGS) -c $*.cpp -o $*.o

costs.o: costs.cpp costs.h
	$(CC) $(CFLAGS) -c $*.cpp -o $*.o

utils.o: utils.cpp utils.h
	$(CC) $(CFLAGS) -c $*.cpp -o $*.o
