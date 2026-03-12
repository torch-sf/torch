#       Makefile for the VETTAM RHD Photoionization Module

RadTrans += rt_ionMomentum.o \
		rt_ionDustRep.o \
		rt_ionisemodule.o \
		rt_ionise.o	\
		rt_ionisedata.o

rt_fillMatrix += rt_ionisemodule.o

rt_source_function_block += rt_ionisemodule.o

rt_ionise.o : rt_ionMomentum.o \
			  rt_ionDustRep.o \
			  rt_ionisemodule.o \
			  rt_ionisedata.o
