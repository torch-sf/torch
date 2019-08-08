#	Makefile for unit RadTrans/RadTransMain/DORIC
RadTrans += rt_data.o rt_init.o\
			calc_ionization.o RadTrans_combine3.o get_ionization.o

rt_init.o: Simulation_data.o
RadTrans.o: calc_ionization.o
RadTrans_combine3.o: calc_ionization.o get_ionization.o
