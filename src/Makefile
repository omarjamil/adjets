##########################################################
# Makefile
# Author: Omar Jamil
SOURCES.cpp = cell.cpp compton.cpp distribution.cpp electrons.cpp numerical.cpp output.cpp paraccess.cpp parameters.cpp photons.cpp run.cpp simulation.cpp synchrotron.cpp volume.cpp 
HEADERS = cell.hpp compton.hpp distribution.hpp electrons.hpp output.hpp paraccess.hpp parameters.hpp photons.hpp simulation.hpp synchrotron.hpp volume.hpp hexa.hpp libs/numerical.hpp libs/physcon.hpp libs/funcobj.hpp
INCLUDES = 
CPP = g++
CFLAGS = 
CPPFLAGS = -O3 -pedantic -std=c++0x 
LIBS = -lstdc++

# the executable file
EFILE= ../bin/aleph
# Object file
OBJDIR= objects
OBJECTS= $(SOURCES.cpp:.cpp=.o)

all: $(EFILE)


# linking to create the executable file
$(EFILE): $(OBJECTS)
	@echo "linking...." 
	$(CPP) $(CPPFLAGS) -o $(EFILE) $(OBJECTS) $(LIBS)


# this declares dependence for object files
# so only the necessary ones
# if the dependence files changed, the .o file
# will be recompiled

cell.o: cell.cpp $(HEADERS) 
	$(CPP) $(CPPFLAGS) -c cell.cpp

compton.o: compton.cpp cell.cpp distribution.cpp electrons.cpp photons.cpp $(HEADERS)
	$(CPP) $(CPPFLAGS) -c compton.cpp

distribution.o: distribution.cpp $(HEADERS)
	$(CPP) $(CPPFLAGS) -c distribution.cpp 

electrons.o: electrons.cpp distribution.cpp $(HEADERS)
	$(CPP) $(CPPFLAGS) -c electrons.cpp

numerical.o: libs/numerical.hpp libs/numerical.cpp
	$(CPP) $(CPPFLAGS) -c libs/numerical.cpp

output.o: output.cpp $(HEADERS)
	$(CPP) $(CPPFLAGS) -c output.cpp

paraccess.o: paraccess.cpp paraccess.hpp 
	$(CPP) $(CPPFLAGS) -c paraccess.cpp

parameters.o: parameters.hpp parameters.cpp 
	$(CPP) $(CPPFLAGS) -c parameters.cpp

photons.o: photons.cpp distribution.cpp $(HEADERS)
	$(CPP) $(CPPFLAGS) -c photons.cpp

run.o: run.cpp
	$(CPP) $(CPPFLAGS) -c run.cpp

simulation.o: simulation.cpp electrons.cpp photons.cpp volume.cpp synchrotron.cpp distribution.cpp $(HEADERS)
	$(CPP) $(CPPFLAGS) -c simulation.cpp

synchrotron.o: synchrotron.cpp electrons.cpp photons.cpp $(HEADERS)
	$(CPP) $(CPPFLAGS) -c synchrotron.cpp

volume.o: volume.cpp synchrotron.cpp photons.cpp electrons.cpp $(HEADERS)
	$(CPP) $(CPPFLAGS) -c volume.cpp




##########################################################
# run "make debug"
# "./debug" will produce "gmon.out"
# "gprof debug" can be used to analyze the output
debug:
	g++ -pg $(OBJECTS) $(LIBS) -o ../$@


# remove the object and executable files
clean:
	rm -f $(EFILE) $(OBJECTS)

spotless:
	rm -f $(EFILE) $(OBJECTS) *~

#
##########################################################