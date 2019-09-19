#include "io_h5_report_xfer_method.h"

int io_h5_report_xfer_method(const int myPE, const hid_t hXferList,
			     const char datasetName[])
{

#ifdef H5_HAVE_PARALLEL
  size_t i, bufspace;
  herr_t err;

/* HDF5 version >= 1.8.8 allows us to report the transfer method
   that was actually used */
# if (H5_VERSION_GE(1,8,8))
  H5D_mpio_actual_io_mode_t mode;
  char mode_str[BUFSIZE];
# endif

/* HDF5 version >= 1.8.10 allows us to report why we did not transfer
   the data using collective I/O */
# if (H5_VERSION_GE(1,8,10))
  char cause_str[BUFSIZE];
  uint32_t local_cause, global_cause;
# endif


# if (H5_VERSION_GE(1,8,8))
  err = H5Pget_mpio_actual_io_mode(hXferList, &mode);
  assert (err >= 0);

  /* There can only be one mode */
  if (myPE == 0) {
    switch (mode) {
    case H5D_MPIO_NO_COLLECTIVE:
      strncpy(mode_str, "H5D_MPIO_NO_COLLECTIVE", BUFSIZE);
      break;
    case H5D_MPIO_CHUNK_INDEPENDENT:
      strncpy(mode_str, "H5D_MPIO_CHUNK_INDEPENDENT", BUFSIZE);
      break;
    case H5D_MPIO_CHUNK_COLLECTIVE:
      strncpy(mode_str, "H5D_MPIO_CHUNK_COLLECTIVE", BUFSIZE);
      break;
    case H5D_MPIO_CHUNK_MIXED:
      strncpy(mode_str, "H5D_MPIO_CHUNK_MIXED", BUFSIZE);
      break;
    case H5D_MPIO_CONTIGUOUS_COLLECTIVE:
      strncpy(mode_str, "H5D_MPIO_CONTIGUOUS_COLLECTIVE", BUFSIZE);
      break;
    default:
      strncpy(mode_str, "Unknown transfer mode", BUFSIZE);
    }
    mode_str[BUFSIZE-1] = '\0';
    printf(" [%s]: Dataset '%s' - transferred in mode %s\n",
	   __FILE__, datasetName, mode_str);
  }
# endif


# if (H5_VERSION_GE(1,8,10))
  err = H5Pget_mpio_no_collective_cause(hXferList, &local_cause, &global_cause);
  assert (err >= 0);

  /* There can be multiple causes */
  if (myPE == 0) {
    if (global_cause != H5D_MPIO_COLLECTIVE) {

      cause_str[0] = '\0';
      bufspace = BUFSIZE;

      /* This macro function appends a cause to cause_str and then reduces the
	 available buffer space by an appropriate amount */
#define ADD_CAUSE_TO_STR(cause)						\
      do {								\
	if (bufspace > 1) {						\
	  strncat(cause_str, cause, bufspace-1);			\
	  bufspace = BUFSIZE - strlen(cause_str);			\
	}								\
      } while(0)


      if (global_cause & H5D_MPIO_SET_INDEPENDENT) {
	ADD_CAUSE_TO_STR("H5D_MPIO_SET_INDEPENDENT ");
      }
      if (global_cause & H5D_MPIO_DATATYPE_CONVERSION) {
	ADD_CAUSE_TO_STR("H5D_MPIO_DATATYPE_CONVERSION ");
      }
      if (global_cause & H5D_MPIO_DATA_TRANSFORMS) {
	ADD_CAUSE_TO_STR("H5D_MPIO_DATA_TRANSFORMS ");
      }
#if (H5_VERSION_LE(1,8,12))
      if (global_cause & H5D_MPIO_SET_MPIPOSIX) {
	ADD_CAUSE_TO_STR("H5D_MPIO_SET_MPIPOSIX ");
      }
#endif
      if (global_cause & H5D_MPIO_NOT_SIMPLE_OR_SCALAR_DATASPACES) {
	ADD_CAUSE_TO_STR("H5D_MPIO_NOT_SIMPLE_OR_SCALAR_DATASPACES ");
      }
#if (H5_VERSION_LE(1,8,12))     
      if (global_cause & H5D_MPIO_POINT_SELECTIONS) {
	ADD_CAUSE_TO_STR("H5D_MPIO_POINT_SELECTIONS ");
      }
#endif
      if (global_cause & H5D_MPIO_NOT_CONTIGUOUS_OR_CHUNKED_DATASET) {
	ADD_CAUSE_TO_STR("H5D_MPIO_NOT_CONTIGUOUS_OR_CHUNKED_DATASET ");
      }
#if (H5_VERSION_LE(1,10,0))
      /* In HDF5 versions starting at 1.10.2(?) there is reportedly
	 H5D_MPIO_NO_COLLECTIVE_MAX_CAUSE instead, not sure whether
	 that has the same meaning or functionality... */
      if (global_cause & H5D_MPIO_FILTERS) {
	ADD_CAUSE_TO_STR("H5D_MPIO_FILTERS ");
      }
#endif

      printf(" [%s]: Dataset '%s' - no collective I/O because %s\n",
	     __FILE__, datasetName, cause_str);
    }
  }
# endif
#endif

  return 0;
}


#ifdef DEBUG_IO

int MPI_File_write_all(MPI_File fh, void *buf, int count,
                       MPI_Datatype datatype, MPI_Status *status)
{
  int myPE, err;
  err = MPI_Comm_rank(MPI_COMM_WORLD, &myPE);
  if (myPE == 0) {
    printf( " [%s]: MPI_File_write_all (collective)\n", __FILE__);
  }
  err = PMPI_File_write_all(fh, buf, count, datatype, status);
  return err;
}

int MPI_File_write_at_all(MPI_File fh, MPI_Offset offset, void *buf,
                          int count, MPI_Datatype datatype,
                          MPI_Status *status)
{
  int myPE, err;
  err = MPI_Comm_rank(MPI_COMM_WORLD, &myPE);
  if (myPE == 0) {
    printf( " [%s]: MPI_File_write_at_all (collective)\n", __FILE__);
  }
  err = PMPI_File_write_at_all(fh, offset, buf, count, datatype, status);
  return err;
}

int MPI_File_write(MPI_File fh, void *buf, int count,
                   MPI_Datatype datatype, MPI_Status *status)
{
  int myPE, err;
  err = MPI_Comm_rank(MPI_COMM_WORLD, &myPE);
  if (myPE == 0) {
    printf( " [%s]: MPI_File_write (independent)\n", __FILE__);
  }
  err = PMPI_File_write(fh, buf, count, datatype, status);
  return err;
}


int MPI_File_write_at(MPI_File fh, MPI_Offset offset, void *buf,
                      int count, MPI_Datatype datatype,
                      MPI_Status *status)
{
  int myPE, err;
  err = MPI_Comm_rank(MPI_COMM_WORLD, &myPE);
  if (myPE == 0) {
    printf( " [%s]: MPI_File_write_at (independent)\n", __FILE__);
  }
  err = PMPI_File_write_at(fh, offset, buf, count, datatype, status);
}

#endif
