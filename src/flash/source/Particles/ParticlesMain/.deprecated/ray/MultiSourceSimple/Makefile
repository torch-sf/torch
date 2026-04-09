Particles += pt_splitRayPoint.o \
		 HEALPixModule.o \
		 pt_blockRayTracePoint.o \
		 pt_initRays.o \
		 pt_finalizeRays.o \
		 pt_generateRaysPoint.o \
		 pt_rayAsyncComm.o \
		 pt_rayAsyncCounter.o

pt_splitRayPoint.o : HEALPixModule.o Particles_rayData.o
HEALPixModule.o : Particles_rayData.o
pt_blockRayTracePoint.o : pt_splitRayPoint.o Particles_rayData.o
pt_rayAsyncCounter.o : pt_rayAsyncComm.o
pt_rayAsyncComm.o : Grid_interface.o
pt_generateRaysPoint.o: rndMT.o
