#	Makefile for the stat+gauss sub-module
SN += SN_data.o SB.o SB_init.o SB_data.o sn_createSNtracer.o rndMT_SN.o SN_init.o
SN.init.o : SN_data.o SB_init.o rndMT_SN.o
SN.o : SN_data.o SB.o SB_data.o
SB.o : SB_data.o rndMT_SN.o
SB_init.o : SB_data.o
