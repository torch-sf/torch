MPIFC ?= mpif90
FC      = $(MPIFC)

LDFLAGS  += -lm $(MUSE_LD_FLAGS)

CODELIB = libflash.a

CODEOBJS = test.o

AR = ar ruv
RANLIB = ranlib
RM = rm

all: $(CODELIB) 

clean:
	$(RM) -f *.o *.a

$(CODELIB): $(CODEOBJS)
	$(RM) -f $@
	$(AR) $@ $(CODEOBJS)
	$(RANLIB) $@

%.o: %.f90
	$(FC) $(FCFLAGS) -c -o $@ $<

