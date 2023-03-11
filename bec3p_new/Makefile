all:	bec3pc

bec3pc:	bec3p.c++ parameters3.h phi-interpolator.c++ phi-interpolator.h
	c++ -O3 bec3p.c++ phi-interpolator.c++ -o bec3pc

phi-test: phi-interpolator.c++ phi-interpolator.h parameters3.h phi-test.c++
	c++ -O3 phi-interpolator.c++ phi-test.c++ -o phi-test

newbec3pc: bec3p.c++ parameters3.h
	c++ -O3 -std=c++11 bec3p.c++ -o bec3pc

phitest: Test_Phi.c++ parameters3.h
	c++ -O3 Test_Phi.c++ -o Test_Phi

readtest: readFile.c++ parameters3.h
	c++ -O3 -std=c++11 readFile.c++ -o readFile
