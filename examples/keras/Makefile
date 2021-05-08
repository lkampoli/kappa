 
 KAPPA_DIR = ../..
 KAPPA_INSTALL = ../../install

 CXX = g++ -lm #-std=c++11
 #CXX = g++ -g -pg -lm -std=c++11 -fopenmp
 #CXX = icpc -std=c++11 -fopenmp
 
 CXX_FLAGS = -I $(KAPPA_INSTALL)/include/kappa++ \
             -I /usr/include 					 \
			 -I /usr/local/include # fdeep

 CXX_LIBS  = -L $(KAPPA_INSTALL)/lib -lkappa++ \
             -L /usr/lib -larmadillo           \
             -L /usr/lib -lyaml-cpp

mixture-sts-shear-keras: mixture-sts-shear-keras.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

%.o : %.cpp
	$(CXX) -c $(CXX_FLAGS)  $<

clean :
	$(RM) *.o $(EXAMPLES) *.out gmon.out *.log

runall:
	./$(EXAMPLES)
