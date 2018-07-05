 
 KAPPA_DIR=../..
 KAPPA_INSTALL=../../install

 CXX       = g++  -g -pg -lm -std=c++11 -fopenmp
 
 CXX_FLAGS = 	-I $(KAPPA_INSTALL)/include/kappa++ 				\
		-I $(KAPPA_DIR)/kappa_deps/x86_64/include			\
		-I $(KAPPA_DIR)/kappa_deps/x86_64/include/armadillo_bits	\
		-I $(KAPPA_DIR)/kappa_deps/x86_64/include/yaml-cpp	

 CXX_LIBS  = 	$(KAPPA_INSTALL)/lib/libkappa++.so				\
		-L $(KAPPA_DIR)/kappa_deps/x86_64/lib -lyaml-cpp -larmadillo	

EXAMPLES  = \
	caching-sts 						\
	heatcond-sts						\
	heatcond-sts-nonrigid					\
	ionized_mixture_creation				\
	mixture_string_creation					\
	mixture-sts-basic					\
	mixture-sts-bulk-nonrigid				\
	mixture-sts-bulk					\
	mixture-sts-bulk-rigid					\
	mixture-sts-bulk_air5					\
	mixture-sts-bulk_REAL_DISTR				\
	mixture-sts-shear					\
	mixture-sts-shear-oxygen				\
	thermodiffusion						\
	mixture-sts-thermo-diffusion				\
	mixture-sts-basic_mixtureName 				\
	mixture_loading						\
	mixture_loading_opt					\
	transport_coefficients					\
	transport_coefficients_air5				\
	mixture-sts-basic_boltzman				\
	TestDiffusion						\
	thermal_conductivity					\
	mixture-sts-thermal_conductivity			\
	mixture-sts-multi-component-diffusion			\
	mixture-sts-multi-component-diffusion_REAL_DISTR	\
	mixture-sts-multi-component-diffusion_air5		\
	mixture-sts-multi-component-diffusion_air5_REAL

all : $(EXAMPLES)

caching-sts : caching-sts.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

heatcond-sts : heatcond-sts.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

heatcond-sts-nonrigid : heatcond-sts-nonrigid.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

ionized_mixture_creation : ionized_mixture_creation.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

mixture_string_creation : mixture_string_creation.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

mixture-sts-basic : mixture-sts-basic.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

mixture-sts-bulk-nonrigid : mixture-sts-bulk-nonrigid.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

mixture-sts-bulk-rigid : mixture-sts-bulk-rigid.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

mixture-sts-bulk : mixture-sts-bulk.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

mixture-sts-bulk_air5 : mixture-sts-bulk_air5.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

mixture-sts-bulk_REAL_DISTR : mixture-sts-bulk_REAL_DISTR.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

mixture-sts-shear : mixture-sts-shear.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

mixture-sts-shear-oxygen : mixture-sts-shear-oxygen.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

thermodiffusion : thermodiffusion.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

mixture-sts-thermo-diffusion : mixture-sts-thermo-diffusion.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

mixture-sts-basic_mixtureName : mixture-sts-basic_mixtureName.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

mixture_loading : mixture_loading.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

mixture_loading_opt : mixture_loading_opt.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

transport_coefficients : transport_coefficients.o
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

mixture-sts-multi-component-diffusion : mixture-sts-multi-component-diffusion.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

mixture-sts-multi-component-diffusion_REAL_DISTR : mixture-sts-multi-component-diffusion_REAL_DISTR.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

mixture-sts-multi-component-diffusion_air5 : mixture-sts-multi-component-diffusion_air5.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

mixture-sts-multi-component-diffusion_air5_REAL : mixture-sts-multi-component-diffusion_air5_REAL.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

%.o : %.cpp
	$(CXX) -c $(CXX_FLAGS)  $<

clean :
	$(RM) *.o $(EXAMPLES) *.out *.txt gmon.out

runall:
	./$(EXAMPLES)
