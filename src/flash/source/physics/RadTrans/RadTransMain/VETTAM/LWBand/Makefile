## Makefile for LW band
RadTrans += rt_lw.o \
            rt_lwdata.o \
            rt_lwmodule.o \
            RayTrace_LW.o

rt_lw.o : rt_ionisedata.o

RayTrace_LW.o : rt_lwdata.o \
                calc_local_block_contributions_sink.o \
                create_cut_block_list_3DRT.o