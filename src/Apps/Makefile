SHELL    = /bin/sh
NAME     = all
MAKEFILE = Makefile

include $(GENIE)/src/make/Make.include
include $(NUPROPEARTH)/src/make/Make.include


BIN_PATH = $(NUPROPEARTH)/bin

TGT = ComputeAttenuation VertexGenerator LeptonPropagator
#TGT = LeptonPropagator

all: $(TGT)

ComputeAttenuation: FORCE 
	$(CXX) $(CXXFLAGS) -c ComputeAttenuation.cxx $(CPP_INCLUDES)
	$(LD) $(LDFLAGS) ComputeAttenuation.o $(LIBRARIES) -o $(BIN_PATH)/ComputeAttenuation

VertexGenerator: FORCE 
	$(CXX) $(CXXFLAGS) -c VertexGenerator.cxx $(CPP_INCLUDES)
	$(LD) $(LDFLAGS) VertexGenerator.o $(LIBRARIES) -o $(BIN_PATH)/VertexGenerator

LeptonPropagator: FORCE 
	$(CXX) $(CXXFLAGS) -c LeptonPropagator.cxx $(CPP_INCLUDES)
	$(LD) $(LDFLAGS) LeptonPropagator.o $(LIBRARIES) -o $(BIN_PATH)/LeptonPropagator

purge:
	$(RM) *.o *~ core

.PHONY: clean
clean:
	@echo "Clean in Apps"
	$(RM) *.o *~ core $(BIN_PATH)/ComputeAttenuation
	$(RM) *.o *~ core $(BIN_PATH)/VertexGenerator
	$(RM) *.o *~ core $(BIN_PATH)/LeptonPropagator

FORCE: