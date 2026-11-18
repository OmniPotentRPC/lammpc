/**
 * Post-step hard-boundary decision. No region geometry lives here: the
 * host callback accepts, rejects, or reflects. A NULL callback is the
 * unset path and does not touch the coordinate buffers.
 */
#include "lammpc_boundary.h"

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

int lammpc_hard_boundary_apply(LammpcHardBoundaryFn fn, void *user,
                               int n_atoms, double *positions,
                               double *velocities,
                               const double *previous_positions,
                               const double *previous_velocities,
                               const LammpcBoundaryBox *box) {
  double *saved_positions;
  double *saved_velocities;
  size_t bytes;
  int decision;

  if (!fn)
    return LAMMPC_BOUNDARY_ACCEPT;
  if (n_atoms < 0 || !box)
    return -1;
  if (n_atoms == 0) {
    decision = fn(user, 0, positions, velocities, box);
    if (decision == LAMMPC_BOUNDARY_ACCEPT ||
        decision == LAMMPC_BOUNDARY_REJECT ||
        decision == LAMMPC_BOUNDARY_REFLECT)
      return decision;
    return -1;
  }
  if (!positions || !velocities || !previous_positions || !previous_velocities)
    return -1;
  if ((size_t)n_atoms > (SIZE_MAX / 3u) / sizeof(double))
    return -1;

  bytes = (size_t)n_atoms * 3u * sizeof(double);
  saved_positions = (double *)malloc(bytes);
  saved_velocities = (double *)malloc(bytes);
  if (!saved_positions || !saved_velocities) {
    free(saved_positions);
    free(saved_velocities);
    return -1;
  }
  memcpy(saved_positions, positions, bytes);
  memcpy(saved_velocities, velocities, bytes);

  decision = fn(user, n_atoms, positions, velocities, box);
  if (decision == LAMMPC_BOUNDARY_ACCEPT) {
    memcpy(positions, saved_positions, bytes);
    memcpy(velocities, saved_velocities, bytes);
  } else if (decision != LAMMPC_BOUNDARY_REFLECT) {
    memcpy(positions, previous_positions, bytes);
    memcpy(velocities, previous_velocities, bytes);
    if (decision != LAMMPC_BOUNDARY_REJECT)
      decision = -1;
  }

  free(saved_positions);
  free(saved_velocities);
  return decision;
}
