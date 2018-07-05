 
 KAPPA_DIR=../..
 KAPPA_INSTALL=../../install

 CXX       = g++ -O3 -pg -lm -std=c++11 -fopenmp
 
 CXX_FLAGS = 	-I $(KAPPA_INSTALL)/include/kappa++ 				\
		-I $(KAPPA_DIR)/kappa_deps/x86_64/include			\
		-I $(KAPPA_DIR)/kappa_deps/x86_64/include/armadillo_bits	\
		-I $(KAPPA_DIR)/kappa_deps/x86_64/include/yaml-cpp	

 CXX_LIBS  = 	-L $(KAPPA_INSTALL)/lib -lkappa++				\
		-L $(KAPPA_DIR)/kappa_deps/x86_64/lib -lyaml-cpp -lopenblas -larmadillo	

EXAMPLES  = basicspeedtest

all : $(EXAMPLES)

basicspeedtest : basicspeedtest.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

%.o : %.cpp
	$(CXX) -c $(CXX_FLAGS)  $<

clean :
	$(RM) *.o $(EXAMPLES) gmon.out
