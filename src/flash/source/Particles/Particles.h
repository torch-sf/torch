!!****h* source/Particles/Particles.h
!!
!! NAME
!! 
!!  Particles.h
!!
!! DESCRIPION
!!
!!  This is the header file for the Particles module.  It includes general
!!  parameters that any routines calling these interfaces needs to
!!  understand.  Interfaces are contained in Particles_interface.F90
!!
!!***

#if 0
 To carry around information about various different particles types,
 the following constants are define to create a data structure 
 particles_typeInfo(PART_TYPE_INFO_SIZE,NPART_TYPES)
 The fields are :
 PART_TYPE - type of particle, PART_MAPMETHOD - mapping method to be
 used in mapping to and from mesh, PART_INTEGRATE_METHOD - time integration
 algorithm, PART_TYPE_BEGIN - the starting point of this particle type, 
 PART_LOCAL - the local number of this type, and 
 PART_LOCAL_MAX - the maximum of this type of particle allowed per proc.
#endif

#define PART_TYPE             1
#define PART_MAPMETHOD        2
#define PART_TYPE_BEGIN       3
#define PART_INITMETHOD       4
#define PART_LOCAL            5
#define PART_LOCAL_MAX        6
#define PART_ADVMETHOD        7
#define PART_TYPE_INFO_SIZE   7

#if 0
   The next few constants relate to the physical quantities associated
   with the particles. Particles are defined with certain attributes.
   One also needs to know the corresponding Grid variable. These constants
   ensure that when a tuple is used with two fields, the first field is 
   always the particle attribute and the second field is the corrensponding
   grid variable index into UNK etc. PART_DS_IND in short for index into
  particle data structure, and GRID_DS_IND is short for index into grid 
  data structure.  The last constant define
  the size of the particles attribute data structure
#endif

#define PART_DS_IND      1
#define GRID_DS_IND      2
#define PART_ATTR_DS_SIZE 2


#if 0
  These constants define the possible mapping methods to and from mesh
#endif
#define QUADRATIC        1
#define WEIGHTED         2
#define NOMAP            3
#define INZONE           4


#if 0
  These constants define the supported integration methods
#endif

#define PASSIVE         323
#define EULER_MAS       364
#define LEAPFROG        366
#define LEAPFROG_COSMO  367

#define EULER_TRA       398
#define RUNGEKUTTA      395
#define MIDPOINT        385
#define ESTI            381
#define CHARGED         382
#define PT_ADVMETH_NONE 383

#define RUNGEKUTTASTATIC 387
#define RUNGEKUTTATHERMO 388
#define MONTECARLO       389

#if 0
  These constants define the supported position initialization methods
#endif

#define LATTICE    678
#define CELLMASS   679
#define REJECTION  680
#define WITH_DENSITY 685
#define FROMFILE   690
#define CUSTOM     691
#define NPARPERZONE   681
#define NPARPERZONEMC 682
#define DPD        867

#define PT_LOGLEVEL_WARN_OPER 200
#define PT_LOGLEVEL_WARN_DATA 300
#define PT_LOGLEVEL_WARN_USE  400

#define PART_EXPAND 321
#define PART_COLLAPSE 344
