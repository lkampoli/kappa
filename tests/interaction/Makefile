
 KAPPA_DIR=../..
 KAPPA_INSTALL=../../install

 CXX       = g++ -O3 -pg -lm -std=c++11 -fopenmp
 
 CXX_FLAGS = 	-I $(KAPPA_INSTALL)/include/kappa++ 				\
#		-I $(KAPPA_DIR)/kappa_deps/x86_64/include		 	\
#		-I $(KAPPA_DIR)/kappa_deps/x86_64/include/armadillo_bits	\
#		-I $(KAPPA_DIR)/kappa_deps/x86_64/include/yaml-cpp	

 CXX_LIBS  = 	-L $(KAPPA_INSTALL)/lib -lkappa++				\
#		-L $(KAPPA_DIR)/kappa_deps/x86_64/lib -lyaml-cpp -larmadillo	

EXAMPLES  = \
	k_exchArrhTest			\
	omegaBCEtest			\
	omegaTest			\
	vssTest				\
	k_diss_ILTtest			\
	k_dissTest			\
        vss_and_neutral_e		\
	test_FHO			\
	cs_dissTest

all : $(EXAMPLES)

k_diss_ILTtest : k_diss_ILTtest.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

k_dissTest : k_dissTest.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

k_exchArrhTest : k_exchArrhTest.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

omegaBCEtest : omegaBCEtest.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

omegaTest : omegaTest.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

vssTest : vssTest.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

vss_and_neutral_e : vss_and_neutral_e.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

test_FHO : test_FHO.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

cs_dissTest : cs_dissTest.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

%.o : %.cpp
	$(CXX) -c $(CXX_FLAGS)  $<

clean :
	$(RM) *.o $(EXAMPLES) *.txt *.ps *.pdf *.xlsx gmon.out *.log
