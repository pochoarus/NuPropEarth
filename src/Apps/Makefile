SHELL    = /bin/sh
NAME     = all
MAKEFILE = Makefile

include $(GENIE)/src/make/Make.include
include $(NUPROPEARTH)/src/make/Make.include


BIN_PATH = $(NUPROPEARTH)/bin


TGT = ComputeAttenuation

all: $(TGT)

ComputeAttenuation: FORCE 
	$(CXX) $(CXXFLAGS) -c ComputeAttenuation.cxx $(CPP_INCLUDES)
	$(LD) $(LDFLAGS) ComputeAttenuation.o $(LIBRARIES) -o $(BIN_PATH)/ComputeAttenuation

purge:
	$(RM) *.o *~ core

.PHONY: clean
clean:
	@echo "Clean in Apps"
	$(RM) *.o *~ core $(BIN_PATH)/ComputeAttenuation

FORCE: