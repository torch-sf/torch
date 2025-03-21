# Makefile for the Hybrid Characteristics radiative transfer module.

RadTrans += RadTrans_hybridCharModule.o		 \
            rt_prepRayTrace.o                  \
            RayTrace_init.o                     \
            RayTrace_finalize.o						\
            RayTrace_data.o							\
	    interpolate_bezier.o                \
            RadTrans_TreeCommData.o             \
            RadTrans_CommunicateGlobalTreeData.o \
            interpolate_bilinear.o 

RadTrans.o : RadTrans_hybridCharModule.o      

RadTrans_hybridCharModule.o : interpolate_bezier.o

RadTrans_CommunicateGlobalTreeData.o : RadTrans_TreeCommData.o

rt_prepRayTrace.o : RayTrace_init.o

#ifdef RAYTRACE_3DRT
RayTrace_init.o : 	rt_init_raytrace_3drt.o 
#endif
