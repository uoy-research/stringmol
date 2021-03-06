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



#############################################################
################## TOP LEVEL TARGETS       ##################
#############################################################

# Notes on the all/debug thing. There's a lot of duplication here, but this 
# helps to keep things clear. We have to run a `cleano' before all/debug because
# otherwise there'll be .o files that may have been built with a different
# config

#Make everything from scratch:
all: $(OUTPUT_DIR)/stringmol $(OUTPUT_DIR)/smspatial

#TODO: Figure out how to set the debug flag to call SDL-devel (or SDL-dev)
#debug: stringmol web_api smspatial
debug: $(OUTPUT_DIR)/stringmol $(OUTPUT_DIR)/smspatial

web_api: bind step

# removes all the object files. Not automatically done by a call to make,
# need to type "make clean" (without the quotes) into the terminal.	
# or "make all clean"
clean:
	rm -f *.o
	rm -f stringmol
	rm -f stringmol_dbg
	rm -f smspatial
	rm -f smspatial_dbg
	rm -f webapi_*.cgi

#just remove the object files:
cleano:
	rm -f *.o
	rm -f ../debug/*.o
	rm -f ../release/*.o

	

#############################################################
################## C OBJECTS	           ##################
#############################################################

COBJECTS = $(OUTPUT_DIR)/hsort.o \
		$(OUTPUT_DIR)/mathutil.o \
		$(OUTPUT_DIR)/memoryutil.o \
		$(OUTPUT_DIR)/microbial_ga.o \
		$(OUTPUT_DIR)/params.o \
		$(OUTPUT_DIR)/randutil.o \
		$(OUTPUT_DIR)/mt19937-2.o \


$(OUTPUT_DIR)/hsort.o: 		hsort.c hsort.h 
	$(C4C) -c hsort.c -o $(OUTPUT_DIR)/hsort.o

$(OUTPUT_DIR)/mathutil.o: 	mathutil.c mathutil.h 
	$(C4C) -c mathutil.c -o $(OUTPUT_DIR)/mathutil.o

$(OUTPUT_DIR)/memoryutil.o: 	memoryutil.c memoryutil.h 
	$(C4C) -c memoryutil.c -o $(OUTPUT_DIR)/memoryutil.o

$(OUTPUT_DIR)/microbial_ga.o: microbial_ga.c microbial_ga.h
	$(C4C) -c microbial_ga.c -o $(OUTPUT_DIR)/microbial_ga.o

$(OUTPUT_DIR)/params.o:		params.c params.h 
	$(C4C) -c params.c -o $(OUTPUT_DIR)/params.o

$(OUTPUT_DIR)/randutil.o: 	randutil.c randutil.h mt19937-2.h
	$(C4C) -c randutil.c -o $(OUTPUT_DIR)/randutil.o

$(OUTPUT_DIR)/mt19937-2.o: 	mt19937-2.c mt19937-2.h
	$(C4C) -c mt19937-2.c -o $(OUTPUT_DIR)/mt19937-2.o



#############################################################
################## C++ OBJECTS	           ##################
#############################################################


CPPOBJECTS = $(OUTPUT_DIR)/lodepng.o \
		$(OUTPUT_DIR)/agents_base.o \
		$(OUTPUT_DIR)/alignment.o \
		$(OUTPUT_DIR)/instructions.o \
		$(OUTPUT_DIR)/rules.o \
		$(OUTPUT_DIR)/setupSM.o \
		$(OUTPUT_DIR)/SMspp.o \
		$(OUTPUT_DIR)/stringmanip.o \
		$(OUTPUT_DIR)/stringPM.o \
		$(OUTPUT_DIR)/tests.o \
		

###################
### PNG library ###
###################

$(OUTPUT_DIR)/lodepng.o: lodepng.cpp lodepng.h
	$(CC) -c lodepng.cpp -o $(OUTPUT_DIR)/lodepng.o


#########################
### Stringmol library ###
#########################

#stringmol cpp files:

$(OUTPUT_DIR)/agents_base.o: agents_base.cpp agents_base.h \
			memoryutil.h	randutil.h	params.h	rules.h
	$(CC) -c agents_base.cpp -o $(OUTPUT_DIR)/agents_base.o

$(OUTPUT_DIR)/alignment.o: 	alignment.cpp alignment.h randutil.h
	$(CC) -c alignment.cpp -o $(OUTPUT_DIR)/alignment.o

$(OUTPUT_DIR)/instructions.o: instructions.cpp instructions.h \
			memoryutil.h randutil.h alignment.h stringmanip.h
	$(CC) -c instructions.cpp -o $(OUTPUT_DIR)/instructions.o

$(OUTPUT_DIR)/rules.o:	rules.cpp rules.h memoryutil.h
	$(CC) -c rules.cpp -o $(OUTPUT_DIR)/rules.o

$(OUTPUT_DIR)/setupSM.o:	setupSM.cpp setupSM.h \
			randutil.h	alignment.h \
			rules.h agents_base.h SMspp.h stringPM.h lodepng.h
	$(CC) -c setupSM.cpp -o $(OUTPUT_DIR)/setupSM.o

$(OUTPUT_DIR)/SMspp.o:	SMspp.cpp SMspp.h
	$(CC) -c SMspp.cpp -o $(OUTPUT_DIR)/SMspp.o

$(OUTPUT_DIR)/stringmanip.o: 	stringmanip.cpp stringmanip.h
	$(CC) -c stringmanip.cpp -o $(OUTPUT_DIR)/stringmanip.o

$(OUTPUT_DIR)/stringPM.o:	stringPM.cpp stringPM.h \
		memoryutil.h	randutil.h	\
        alignment.h		stringmanip.h	instructions.h	\
		rules.h	agents_base.h	SMspp.h 
	$(CC) -c stringPM.cpp -o $(OUTPUT_DIR)/stringPM.o

$(OUTPUT_DIR)/tests.o:	tests.cpp tests.h
	$(CC) -c tests.cpp -o $(OUTPUT_DIR)/tests.o

	
	
	
	
#############################################################
################## LIBS and EXECUTABLES    ##################
#############################################################


#############################################################
# spatial stringmol

$(OUTPUT_DIR)/smspatial: 	$(OUTPUT_DIR)/smspatial.o \
			$(COBJECTS) $(CPPOBJECTS) 	
	$(CC) -o smspatial -s -lm 	\
	        -funroll-loops -fomit-frame-pointer \
	        $(OUTPUT_DIR)/smspatial.o \
			$(COBJECTS) $(CPPOBJECTS) 
#			-lSDL 
		
		
###########################################################
# stringmol
$(OUTPUT_DIR)/stringmol: 	$(OUTPUT_DIR)/stringmol.o \
			$(COBJECTS) $(CPPOBJECTS)
	$(CC) -o $(OUTPUT_DIR)/stringmol -lm $(OUTPUT_DIR)/stringmol.o \
			$(COBJECTS) $(CPPOBJECTS)\
			
				

###########################################################
# stringmol web api			
			
bind: 	webapi_bind.o webapi_util.o stringPM.o agents_base.o rules.o SMspp.o setupSM.o \
			$(COBJECTS) \
			stringmanip.o	alignment.o	instructions.o tests.o lodepng.o
	$(CC) -o webapi_bind.cgi -lm webapi_bind.o webapi_util.o stringPM.o agents_base.o rules.o SMspp.o setupSM.o \
			$(COBJECTS) \
			stringmanip.o	alignment.o	instructions.o tests.o lodepng.o
			
			
step: 	webapi_step.o webapi_util.o stringPM.o agents_base.o rules.o SMspp.o setupSM.o \
			$(COBJECTS) \
			stringmanip.o	alignment.o	instructions.o tests.o  lodepng.o
	$(CC) -o webapi_step.cgi -lm webapi_step.o webapi_util.o stringPM.o agents_base.o rules.o SMspp.o setupSM.o \
			$(COBJECTS) \
			stringmanip.o	alignment.o	instructions.o tests.o	 lodepng.o		
		


#########################
### stringmol objects ###
#########################


$(OUTPUT_DIR)/stringmol.o:  	stringmol.cpp \
			$(COBJECTS) $(CPPOBJECTS)
	$(CC) -c stringmol.cpp -o  	$(OUTPUT_DIR)/stringmol.o
			

$(OUTPUT_DIR)/smspatial.o:  	smspatial.cpp \
			$(COBJECTS) $(CPPOBJECTS)
	$(CC) -c smspatial.cpp 	-o $(OUTPUT_DIR)/smspatial.o

#########################
### web_api objects ###
#########################


webapi_util.o:	webapi_util.cpp webapi_util.h stringPM.h \
		randutil.h	params.h hsort.h\
        alignment.h	setupSM.h\
		rules.h	agents_base.h	SMspp.h 
	$(CC) -c webapi_util.cpp
	
	
	
webapi_bind.o:  webapi_bind.cpp webapi_util.h\
			randutil.h	alignment.h \
			rules.h agents_base.h SMspp.h stringPM.h \
			tests.h webapi_util.h
	$(CC) -c webapi_bind.cpp 

webapi_step.o:  webapi_step.cpp webapi_util.h\
			randutil.h	alignment.h \
			rules.h agents_base.h SMspp.h stringPM.h \
			tests.h webapi_util.h
	$(CC) -c webapi_step.cpp 		

