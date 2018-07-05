
 KAPPA_DIR=../..
 KAPPA_INSTALL=../../install

 CXX       = g++ -O3 -pg -lm -std=c++11 -fopenmp
 
 CXX_FLAGS = 	-I $(KAPPA_INSTALL)/include/kappa++ 				\
		-I $(KAPPA_DIR)/kappa_deps/x86_64/include		 	\
		-I $(KAPPA_DIR)/kappa_deps/x86_64/include/armadillo_bits	\
		-I $(KAPPA_DIR)/kappa_deps/x86_64/include/yaml-cpp	

 CXX_LIBS  = 	-L $(KAPPA_INSTALL)/lib -lkappa++				\
		-L $(KAPPA_DIR)/kappa_deps/x86_64/lib -lyaml-cpp -larmadillo	

EXAMPLES  = \
	cs_dissTest 			\
	H2H_omegatest			\
	interaction_speedtest		\
	interactionTest			\
	k_exchArrhTest			\
	k_exchTest			\
	k_VTTest			\
	k_VTTest_LC			\
	omegaBCEtest			\
	omega_nc_test			\
	omegaTest			\
	vssTest				\
	k_diss_ILTtest			\
	k_dissTest			\
	k_dissTest_LC			\
        vss_and_neutral_e		\
	R_VT				\
	R_VV				\
	R_VV_prime			\
	R_diss-rec			\
	R_exch				\
	probability_transition_VT	\
	test_FHO

all : $(EXAMPLES)

cs_dissTest : cs_dissTest.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

H2H_omegatest : H2H_omegatest.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

interaction_speedtest : interaction_speedtest.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

interactionTest : interactionTest.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

k_diss_ILTtest : k_diss_ILTtest.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

k_dissTest : k_dissTest.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

k_dissTest_LC : k_dissTest_LC.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

k_exchArrhTest : k_exchArrhTest.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

k_exchTest : k_exchTest.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

k_VTTest : k_VTTest.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

k_VTTest_LC : k_VTTest_LC.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

omegaBCEtest : omegaBCEtest.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

omega_nc_test : omega_nc_test.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

omegaTest : omegaTest.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

vssTest : vssTest.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

vss_and_neutral_e : vss_and_neutral_e.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

R_VT : R_VT.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

R_VV : R_VV.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

R_VV_prime : R_VV_prime.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

R_diss-rec : R_diss-rec.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

R_exch : R_exch.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

#probability_transition_VT : probability_transition_VT.o
#	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

test_FHO : test_FHO.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

%.o : %.cpp
	$(CXX) -c $(CXX_FLAGS)  $<

clean :
	$(RM) *.o $(EXAMPLES) *.txt *.ps *.pdf *.xlsx gmon.out

