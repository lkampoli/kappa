 
 KAPPA_DIR = ../..
 KAPPA_INSTALL = ../../install

 CXX = g++ -g -pg -lm -std=c++11 -fopenmp
 
 CXX_FLAGS = -I $(KAPPA_INSTALL)/include/kappa++             \
#		-I $(KAPPA_DIR)/kappa_deps/x86_64/include                \
#		-I $(KAPPA_DIR)/kappa_deps/x86_64/include/armadillo_bits \
#		-I $(KAPPA_DIR)/kappa_deps/x86_64/include/yaml-cpp

 CXX_LIBS  = 	$(KAPPA_INSTALL)/lib/libkappa++.so \
#		-L $(KAPPA_DIR)/kappa_deps/x86_64/lib -lyaml-cpp -larmadillo -lopenblas -llapack

EXAMPLES  = 											 \
	mixture_string_creation          \
	mixture-sts-basic                \
	mixture-sts-bulk_air5            \
	mixture-sts-shear                \
	mixture-sts-bulk                 \
	mixture-sts-shear-bulk-thermal   \
	thermodiffusion                  \
	mixture-sts-thermo-diffusion     \
	transport_coefficients_air5      \
	mixture-sts-basic_boltzman       \
	TestDiffusion                    \
	thermal_conductivity             \
	mixture-sts-thermal_conductivity 
#	ionized_mixture_creation         \


all : $(EXAMPLES)

ionized_mixture_creation : ionized_mixture_creation.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

mixture_string_creation : mixture_string_creation.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

mixture-sts-basic : mixture-sts-basic.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

mixture-sts-bulk_air5 : mixture-sts-bulk_air5.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

mixture-sts-shear : mixture-sts-shear.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

mixture-sts-bulk : mixture-sts-bulk.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

mixture-sts-shear-bulk-thermal : mixture-sts-shear-bulk-thermal.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

thermodiffusion : thermodiffusion.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

mixture-sts-thermo-diffusion : mixture-sts-thermo-diffusion.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

transport_coefficients_air5 : transport_coefficients_air5.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

mixture-sts-basic_boltzman : mixture-sts-basic_boltzman.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

TestDiffusion : TestDiffusion.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

thermal_conductivity : thermal_conductivity.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

mixture-sts-thermal_conductivity : mixture-sts-thermal_conductivity.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

%.o : %.cpp
	$(CXX) -c $(CXX_FLAGS)  $<

clean :
	$(RM) *.o $(EXAMPLES) *.out *.txt gmon.out *.log

runall:
	./$(EXAMPLES)
