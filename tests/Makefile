#############################################################
################## VARIABLES               ##################
#############################################################


# which compiler to use 
C4C = gcc -O3 -Wall -Wunused
CC = g++ -O3 -Wall -Wunused -std=c++11

#TODO: Handle profiling better...
debug: C4C = gcc -O0 -Wall -Wunused -DDEBUG -g 
debug: CC = g++ -O0 -Wall -Wunused -std=c++11 -DDEBUG -g

#If you want profiling, comment out the above and replace with: 
#debug: C4C = gcc -O0 -Wall -Wunused -DDEBUG -g -pg
#debug: CC = g++ -O0 -Wall -Wunused -std=c++11 -DDEBUG -g -pg

#build directories:
DEBUG_DIR =  ../debug
RELEASE_DIR =  ../release

#TODO: There's a lot of repetition of $(OUTPUT_DIR) atm...!!!
OUTPUT_DIR = $(RELEASE_DIR)
debug: OUTPUT_DIR = $(DEBUG_DIR)



RP="../debug"

#TODO: can't use wildcards e.g. ../release/*.o because of multiple 'main's... fix!
tests: 
	g++ -std=gnu++11 -Wall  -DDEBUG -g -o test  *.cpp ${RP}/mt19937-2.o ${RP}/randutil.o \
  	${RP}/SMspp.o ${RP}/stringPM.o ${RP}/agents_base.o ${RP}/rules.o ${RP}/alignment.o \
  	${RP}/params.o ${RP}/memoryutil.o ${RP}/instructions.o ${RP}/stringmanip.o \
  	${RP}/hsort.o




