ROOTCONFIG   := root-config

ARCH         := $(shell $(ROOTCONFIG) --arch)
PLATFORM     := $(shell $(ROOTCONFIG) --platform)


ROOTCFLAGS   := $(shell $(ROOTCONFIG) --cflags)
ROOTLDFLAGS  := $(shell $(ROOTCONFIG) --ldflags)
ROOTLIBS     := $(shell $(ROOTCONFIG) --libs)
ROOTGLIBS    := $(shell $(ROOTCONFIG) --glibs)
HASTHREAD    := $(shell $(ROOTCONFIG) --has-thread)
ROOTINC      := $(shell $(ROOTCONFIG) --incdir)


CC            = gcc
CXX           = g++
LD            = g++
OutPutOpt     = -o # keep whitespace after "-o"

#
ifeq ($(ARCH),linux)
CXXFLAGS      = -O2 -Wall -fPIC
LDFLAGS       = -O2
SOFLAGS       = -shared
endif

#
ifeq ($(ARCH),linuxx8664gcc)
CXXFLAGS      = -O2 -Wall -fPIC
LDFLAGS       = -O2 -Wl,--no-as-needed $(ROOTGLIBS)
SOFLAGS       = -shared
endif

#
ifeq ($(ARCH),win32gcc)
CXXFLAGS      = -O -pipe -Wall -Woverloaded-virtual -I/usr/X11R6/include
LDFLAGS       = -O -Wl,--enable-auto-import -Wl,--enable-runtime-pseudo-reloc \
		-L/usr/X11R6/lib
SOFLAGS       = -shared -D_DLL -Wl,--export-all-symbols
endif


CXXFLAGS     += $(ROOTCFLAGS)
LDFLAGS      += $(ROOTLDFLAGS) 
LIBS          = $(ROOTLIBS) $(SYSLIBS)
GLIBS         = $(ROOTGLIBS) $(SYSLIBS)


#------------------------------------------------------------------------------

pro           = drs_01

OBJS          = $(pro).o  

#------------------------------------------------------------------------------

#CLHEPINC =  $(CLHEP_PARAM_PATH)/include
#CLHEPLIB =  $(CLHEP_PARAM_PATH)/lib

#CLHEPINC :=  $(shell clhep-config --include)
#CLHEPLIB :=  $(shell clhep-config --libs)

#INCLUDES_C =    $(CLHEPINC)  

#LIBS2  = $(CLHEPLIB) -lCLHEP


#------------------------------------------------------------------------------

all:            $(pro)


$(pro):       $(OBJS)
		$(LD)  $(LDFLAGS) $(OBJS) $(LIBS) -lMathMore $(GLIBS)  -o $(pro)



clean:
		@rm -f *.o core.* *~

#------------------------------------------------------------------------------


%.o:%.cc
	$(CXX) $(CXXFLAGS) $(INCLUDES_C:%=%) -c $<

%.o:%.c
	$(CC) $(CXXFLAGS) $(INCLUDES_C:%=%) -c $<




