Heat += CoolVars.o heatCool.o HeatCoolInterface.o

HeatCoolInterface.o: OdeData.o
CoolVars.o: OdeData.o
heatCool.o: CoolVars.o OdeInt.o Cool_data.o find.o
RadHeat.o: heatCool.o
