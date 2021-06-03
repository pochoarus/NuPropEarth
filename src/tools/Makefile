SHELL    = /bin/sh
NAME     = all
MAKEFILE = Makefile

include $(GENIE)/src/make/Make.include
include $(NUPROPEARTH)/src/make/Make.include


BIN_PATH = $(NUPROPEARTH)/bin


TGT = BuildEarth

all: $(TGT)

BuildEarth: FORCE 
	$(CXX) $(CXXFLAGS) -c BuildEarth.cxx $(CPP_INCLUDES)
	$(LD) $(LDFLAGS) BuildEarth.o $(LIBRARIES) -o $(BIN_PATH)/BuildEarth

purge:
	$(RM) *.o *~ core

.PHONY: clean
clean:
	@echo "Clean in Apps"
	$(RM) *.o *~ core $(BIN_PATH)/BuildEarth

FORCE: