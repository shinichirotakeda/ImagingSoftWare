TARGET        = sbp2d_norm
OBJECTS = SBPComptonImaging.o TwoHitComptonEvent.o OneEventResponse.o

##############################################################################

#------------------------------------------------------------#

#CXX           = g++
CXX           = g++ -std=c++11
CC            = gcc

ObjSuf        = o
CXXSrcSuf     = cc

OutPutOpt     = -o 

#CXXFLAGS      = -O -Wall -fPIC -fopenmp -I. 
CXXFLAGS      = -O -Wall -fPIC -I. 

LD            = g++ -bind_at_load
LDFLAGS       = -O

#------------------------------------------------------------#

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)

NEXT_CLI_INC = $(HOME)/include/anl
NEXT_CLI_LNK = $(HOME)/lib/libANLNext.dylib

BOOST_INC_DIR = /usr/local/include
BOOST_LIB_DIR = /usr/local/lib
BOOST_LIB = $(BOOST_LIB_DIR)/libboost_system-mt.dylib $(BOOST_LIB_DIR)/libboost_filesystem-mt.dylib

SYSLIB = -lreadline

CXXFLAGS += -I$(NEXT_CLI_INC) -I$(BOOST_INC_DIR) $(ROOTCFLAGS)

#LIBS = $(COM_CLI_LNK) $(ROOTLIBS) $(ROOTGLIBS) $(SYSLIB) -lgomp
LIBS = $(NEXT_CLI_LNK) $(ROOTLIBS) $(ROOTGLIBS) $(SYSLIB) $(BOOST_LIB)



##############################################################################

all: $(TARGET)

.cc.o:
	$(CXX) -O2 $(CXXFLAGS) -c $< -I. \
	

$(TARGET) : %: %.o $(OBJECTS)
	$(CXX) -bind_at_load -o $@ $^ \
	$(LIBS)

clean:
	@rm -f $(TARGET) *.o *~ 

