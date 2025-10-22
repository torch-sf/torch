# Makefile for the Hybrid Characteristics radiative transfer module.

RadTrans += RadTrans_hybridCharModule.o		 \
            rt_prepRayTrace.o                  \
            raytrace_data.o							\
	    interpolate_bezier.o                \
            RadTrans_TreeCommData.o             \
            RadTrans_CommunicateGlobalTreeData.o \
            interpolate_bilinear.o 

RadTrans.o : RadTrans_hybridCharModule.o      

RadTrans_hybridCharModule.o : interpolate_bezier.o

RadTrans_CommunicateGlobalTreeData.o : RadTrans_TreeCommData.o


#ifdef RAYTRACE_3DRT
raytrace_init.o : 	rt_init_raytrace_3drt.o 
#endif
