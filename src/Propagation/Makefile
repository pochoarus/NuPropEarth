SHELL    = /bin/sh
NAME     = all
MAKEFILE = Makefile

include $(GENIE)/src/make/Make.include
include $(NUPROPEARTH)/src/make/Make.include

PACKAGE           = Propagation
DICTIONARY        = _ROOT_DICT_Propagation
LIBNAME           = libGPropagation

all     : rootcint lib

include $(NUPROPEARTH)/src/make/Make.std-package-targets

FORCE: