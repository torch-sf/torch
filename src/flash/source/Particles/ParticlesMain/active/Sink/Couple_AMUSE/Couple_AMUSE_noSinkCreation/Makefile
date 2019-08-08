Particles += dens_removal.o sampleIMF.o \
             gaussian.o sort.o \
             Particles_sinkMakeStars.o \
             Particles_starForm.o \
             stochastic.o

dens_removal.o: sampleIMF.o gaussian.o qsort2D.o ut_qsort.o
Particles_sinkMakeStars.o: sampleIMF.o
Particles_starForm.o: sampleIMF.o sort.o gaussian.o
