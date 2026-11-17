#ifndef LAMMPC_BOUNDARY_EMBED_H
#define LAMMPC_BOUNDARY_EMBED_H

/**
 * Internal hooks used by the embed backend to register fix lammpc/boundary
 * on a live LAMMPS instance. Not part of the installed C ABI. The fix is
 * created only while a host callback is armed for a run.
 */
#include "lammpc_boundary.h"

#ifdef __cplusplus
extern "C" {
#endif

int lammpc_embed_boundary_attach(void *lmp, LammpcHardBoundaryFn fn,
                                 void *user);
int lammpc_embed_boundary_detach(void *lmp);
void lammpc_embed_boundary_disarm(void *lmp);
const char *lammpc_embed_boundary_error(void);

#ifdef __cplusplus
}
#endif

#endif /* LAMMPC_BOUNDARY_EMBED_H */
