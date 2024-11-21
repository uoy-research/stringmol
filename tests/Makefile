#############################################################
################## VARIABLES               ##################
#############################################################


# which compiler to use 
CC = g++ -O3 -Wall -Wunused -std=c++11

# relative path to debug objects
RP="../debug"

TESTFILES = stringPM_test.cpp \
	    rngTest.cpp \
	    

#TODO: can't use wildcards e.g. ../release/*.o because of multiple 'main's... fix!
tests: 
	g++ -std=gnu++11 -Wall  -DDEBUG -g -o test  *.cpp ${RP}/mt19937-2.o ${RP}/randutil.o \
  	${RP}/SMspp.o ${RP}/stringPM.o ${RP}/agents_base.o ${RP}/rules.o ${RP}/alignment.o \
  	${RP}/params.o ${RP}/memoryutil.o ${RP}/instructions.o ${RP}/stringmanip.o \
  	${RP}/hsort.o ${RP}/opcodes.o\
  	catch.o


catch.o:  	RunTest.cpp
	$(CC) -c RunTest.cpp -o catch.o



