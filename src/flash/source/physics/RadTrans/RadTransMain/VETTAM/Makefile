#	Makefile for the VETTAM RHD Module

RadTrans += RadTrans_data.o \
	rt_petsc.o \
    rt_fillMatrix.o \
    rt_dustTerms.o \
	rt_dustTemperature.o \
	rt_sinkInject.o \
	rt_setOpacity.o \
	rt_setEddTensor.o \
	rt_sinkHydro.o

	
#These dependencies should be invoked by default, but no harm in adding them here.
RadTrans_init.o : RadTrans_data.o
RadTrans.o : RadTrans_data.o
rt_petsc.o : RadTrans_data.o \
	       rt_fillMatrix.o
rt_fillMatrix.o : RadTrans_data.o \
	       rt_setOpacity.o
rt_dustTerms.o : RadTrans_data.o
rt_dustTemperature.o : RadTrans_data.o