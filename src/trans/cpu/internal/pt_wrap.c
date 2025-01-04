#include <mpi.h>
#include "progress_helper.h"

void pt_reqset_register_f(int *count,int *array_of_requests,int *flags,int *gid,int *status)
{
  
  int i;
  MPI_Request *myreqs,*pr;
  myreqs = malloc(*count * sizeof(MPI_Request));
  pr = myreqs;
  for(i=0;i< *count;i++)
    *pr++ = MPI_Request_f2c(array_of_requests[i]);

  pt_reqset_register(*count,myreqs,*flags,gid,status);
}
