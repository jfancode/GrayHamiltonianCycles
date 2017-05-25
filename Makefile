# The following Makefile was created with edits from thc's lecture notes
COMPILER = clang++ -O1
CILKARGS = -fcilkplus -lcilkrts
UTILFLAGS = -lm

GrayHamiltonian: GrayHamiltonian.cpp util.cpp
	$(COMPILER) $(CILKARGS) -o $@ $^ $(UTILFLAGS)

GrayHamiltoniannocilk: GrayHamiltonian.cpp util.cpp
	$(COMPILER) $(CILKARGS) -o $@ -DNOCILK $^ $(UTILFLAGS)

