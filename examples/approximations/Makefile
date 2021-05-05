 
 KAPPA_DIR=../..
 KAPPA_INSTALL=../../install

 CXX       = g++ -O3 -pg -lm -std=c++11 -fopenmp

 CXX_FLAGS = -I $(KAPPA_INSTALL)/include/kappa++ \
	     -I /usr/include

 CXX_LIBS  = -L $(KAPPA_INSTALL)/lib -lkappa++ \
	     -L /usr/lib -larmadillo           \
	     -L /usr/lib -lyaml-cpp

 EXAMPLES  = cvibr_multit \
	     cvibr_sts    \
     	     cvtr_sts     \
     	     cvrot_sts
	
all : $(EXAMPLES)

cvibr_multit : cvibr_multit.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

cvibr_sts : cvibr_sts.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

cvtr_sts : cvtr_sts.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

cvrot_sts : cvrot_sts.o
	$(CXX) -o $@ $(CXX_FLAGS) $^ $(CXX_LIBS)

%.o : %.cpp
	$(CXX) -c $(CXX_FLAGS)  $<

clean :
	$(RM) *.o $(EXAMPLES) *.txt gmon.out *.log
