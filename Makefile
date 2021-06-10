SHELL    = /bin/sh
NAME     = all
MAKEFILE = Makefile

include $(GENIE)/src/make/Make.include

BUILD_TARGETS = Flux Driver Propagation Apps tools

all:     $(BUILD_TARGETS)

Flux: FORCE
	@echo " "
	@echo "** Building Flux..."
	test -d lib || mkdir lib
	cd $(NUPROPEARTH)/src/Flux && $(MAKE)
	cd $(NUPROPEARTH)

Driver: FORCE
	@echo " "
	@echo "** Building Driver..."
	test -d lib || mkdir lib
	cd $(NUPROPEARTH)/src/Driver && $(MAKE)
	cd $(NUPROPEARTH)

Propagation: FORCE
	@echo " "
	@echo "** Building Propagation..."
	test -d lib || mkdir lib
	cd $(NUPROPEARTH)/src/Propagation && $(MAKE)
	cd $(NUPROPEARTH)

Apps: FORCE
	@echo " "
	@echo "** Building Apps..."
	test -d bin || mkdir bin && \
	cd $(NUPROPEARTH)/src/Apps && $(MAKE)
	cd $(NUPROPEARTH)

tools: FORCE
	@echo " "
	@echo "** Building tools..."
	test -d bin || mkdir bin
	cd ${NUPROPEARTH}/src/tools && $(MAKE)
	cd ${NUPROPEARTH}




clean: FORCE
	@echo " "
	@echo "** Cleaning all.. "
	cd $(NUPROPEARTH)/src/Flux && $(MAKE) clean
	cd $(NUPROPEARTH)/src/Driver && $(MAKE) clean
	cd $(NUPROPEARTH)/src/Propagation && $(MAKE) clean
	cd $(NUPROPEARTH)/src/Apps && $(MAKE) clean
	rm -f $(NUPROPEARTH)/lib/* $(NUPROPEARTH)/bin/*
	cd $(NUPROPEARTH)


FORCE:
