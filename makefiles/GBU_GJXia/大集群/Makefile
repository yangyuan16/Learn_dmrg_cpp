# Lines starting with the pound sign are comments.
#
# This is one of two  options you might need to tweek.

EXECUTABLE=runDMRG
#slurm_cpus_per_node=$(abcde)
# You can modify the below as well, but probably won't need to

# CC is for the name of the C++ compiler. CPPFLAGS denotes pre-processor
# flags, such as -l options. CFLAGS denotes flags for the C compiler.
# CXXFLAGS denotes flags for the C++ compiler. You may additional
# setting here, such as PFLAFS, if you are using other languages such
# as Pascal.

CC=icc -diag-disable=10441
CPPFLAGS= -O2 -xHost -I$(GSL_HOME)/include/ -Dslurm_cpus_per_node=$(slurm_cpus_per_node) 
CFLAGS=-g -Wall -O2 -xHost  -I$(GSL_HOME)/include/ -Dslurm_cpus_per_node=$(slurm_cpus_per_node)
CXXFLAGS=$(CFLAGS)
COMPILE=$(CC) $(CPPFLAGS) $(CXXFLAGS) -c

SRCS:=$(wildcard *.cpp)
OBJS:=$(patsubst %.cpp,%.o,$(SRCS))
DEPS:=$(patsubst %.cpp,%.d,$(SRCS))
#LIBS=-lstdc++    -L/home/xinlu/usr/lib -lgsl  -mkl
LIBS=-lstdc++    -L$(GSL_HOME)/lib  -lgsl  -mkl

# "all" is the default target. Simply make it point to myprogram.

all:$(EXECUTABLE)

# Define the components of the program, and how to link them together.
# These components are defined as dependencies; that is,
# they must be made up-to-date before the code is linked.

$(EXECUTABLE): $(DEPS) $(OBJS)
	$(CC) -o $(EXECUTABLE) $(OBJS) $(LIBS)

# Specify that the dependency files depend on the C++ cource files.

%.d:%.cpp
	$(CC) -MM $(CPPFLAGS) $< > $@
	$(CC) -MM $(CPPFLAGS) $< | sed s/\\.o/\\.d/   >> $@

# Specify that all .o files depend on .cpp files, and indicate how
# the .cpp files are converted(compiled) to the .o files. 

%.o:%.cpp
	$(COMPILE) -o $@ $<

clean:
	-rm $(OBJS) $(DEPS) $(EXECUTABLE) 

explain:
	@echo "The following information represents your prgram"
	@echo "Final executable name: $(EXECUTABLE)"
	@echo "Source files: $(SRCS)"
	@echo "Object files: $(OBJS)"
	@echo "Dependency files: $(DEPS)"

depend:$(DEPS)
	@echo "Dependencies are now up-to-date."

-include $(DEPS)

