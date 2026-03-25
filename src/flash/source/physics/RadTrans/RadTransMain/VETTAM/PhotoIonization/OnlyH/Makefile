RadTrans += rt_ionHeatCool.o \
		rt_ionRateEqs.o 

rt_ionise.o : rt_ionHeatCool.o \
			rt_ionRateEqs.o \
			rt_sinkInject.o \
			rt_setOpacity.o

rt_ionRateEqs.o : rt_ionisemodule.o \
			  	  rt_ionisedata.o
