#ifndef PROGRESS_HELPER_HEADER_HAS_BEEN_INCLUDED
#define PROGRESS_HELPER_HEADER_HAS_BEEN_INCLUDED

#include <mpi.h>

typedef enum pt_request_set_flags_e {
    PT_REQSET_NON_PERSISTENT = 1,
    PT_REQSET_STATUSES_IGNORE = 2,
    /* Not to be set by the callee */
    PT_REQSET_ACTIVE = 4,
    PT_REQSET_COMPLETED = 8,
} pt_request_set_flags_e;

void pt_reqset_wait(int *gid, int *status);
void pt_reqset_test(int *gid, int* flag, int *status);
void pt_reqset_start(int *gid,int *status);
int pt_reqset_unregister(int* gid);
void pt_reqset_register(int count, MPI_Request* array_of_requests, int flags, int* gid, int *status);

#endif  /* PROGRESS_HELPER_HEADER_HAS_BEEN_INCLUDED */
