 
 KAPPA_DIR=../..
 KAPPA_INSTALL=../../install

 CXX       = g++ -O3 -pg -lm -std=c++11 -fopenmp
 
 CXX_FLAGS = -I $(KAPPA_INSTALL)/include/kappa++ \
	     -I /usr/include

 CXX_LIBS  = -L $(KAPPA_INSTALL)/lib -lkappa++ \
	     -L /usr/lib -larmadillo           \
             -L /usr/lib -lyaml-cpp

 EXAMPLES  = basicTest

all : $(EXAMPLES)

basicTest : basicTest.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

%.o : %.cpp
	$(CXX) -c $(CXX_FLAGS)  $<

clean :
	$(RM) *.o $(EXAMPLES) gmon.out *.log
