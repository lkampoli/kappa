 
 KAPPA_DIR=../..
 KAPPA_INSTALL=../../install

 CXX       = g++ -O3 -pg -lm -std=c++11 -fopenmp
 
 CXX_FLAGS = 	-I $(KAPPA_INSTALL)/include/kappa++ 				\
		-I $(KAPPA_DIR)/kappa_deps/x86_64/include		 	\
		-I $(KAPPA_DIR)/kappa_deps/x86_64/include/armadillo_bits	\
		-I $(KAPPA_DIR)/kappa_deps/x86_64/include/yaml-cpp	

 CXX_LIBS  = 	$(KAPPA_INSTALL)/lib/libkappa++.so 				\
		-L$(KAPPA_DIR)/kappa_deps/x86_64/lib -lyaml-cpp -larmadillo	

EXAMPLES  = \
	dumpspectrumTest	\
	particleProperties	\
	particleTest		\
	particle_loading		
	

all : $(EXAMPLES)

dumpspectrumTest : dumpspectrumTest.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

particleProperties : particleProperties.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

particleTest : particleTest.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

particle_loading : particle_loading.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

%.o : %.cpp
	$(CXX) -c $(CXX_FLAGS)  $<

clean :
	$(RM) *.o $(EXAMPLES) *.csv *.out gmon.out

