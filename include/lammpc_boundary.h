#ifndef LAMMPC_BOUNDARY_H
#define LAMMPC_BOUNDARY_H

/**
 * @file lammpc_boundary.h
 * @brief Host hard-boundary decision applied after an integration step.
 *
 * lammpc does not implement a region. The host does, and reports one of
 * three decisions. This matches the hard-region rule used by a sampler:
 * a configuration inside the closed region is accepted, a predicate
 * failure restores the previous contained state, and a planar crossing
 * is kept only when the host has already reflected positions and
 * velocities in the buffers passed to the callback.
 *
 * The callback is optional. `lammpc_hard_boundary_apply()` returns
 * immediately when the function pointer is NULL and does not read or
 * write the coordinate buffers.
 */

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Decision returned by a hard-boundary callback. */
typedef enum LammpcBoundaryDecision {
  /** @brief Keep the integrator state from after the step. */
  LAMMPC_BOUNDARY_ACCEPT = 0,
  /** @brief Restore positions and velocities from before the step. */
  LAMMPC_BOUNDARY_REJECT = 1,
  /**
   * @brief Keep positions and velocities the callback wrote in place.
   * Writes are ignored for every other decision.
   */
  LAMMPC_BOUNDARY_REFLECT = 2
} LammpcBoundaryDecision;

/**
 * @brief Simulation box passed to the callback, in LAMMPS metal units.
 *
 * Lengths are angstroms. `lo` and `hi` are the orthogonal edges. `xy`,
 * `xz`, and `yz` are the restricted-triclinic tilt factors (zero when
 * `triclinic` is 0). `periodic` is 1 on each axis LAMMPS treats as periodic.
 */
typedef struct LammpcBoundaryBox {
  double lo[3];
  double hi[3];
  double xy;
  double xz;
  double yz;
  int periodic[3];
  int triclinic;
} LammpcBoundaryBox;

/**
 * @brief Host decision for one post-step configuration.
 *
 * `positions` and `velocities` have `3 * n_atoms` entries in atom-tag
 * order (tag 1 is index 0). Positions are angstroms and velocities are
 * angstroms per picosecond, which is how a metal-unit LAMMPS instance
 * stores them. The host may write those buffers only when it returns
 * `LAMMPC_BOUNDARY_REFLECT`.
 *
 * A planar face is the reflect path: the host mirrors the crossing
 * coordinate and the outward normal component of velocity, including any
 * projection into an active constraint tangent space. A predicate with no
 * penetration depth is the reject path. A point in the closed region is
 * the accept path.
 */
typedef int (*LammpcHardBoundaryFn)(void *user, int n_atoms, double *positions,
                                    double *velocities,
                                    const LammpcBoundaryBox *box);

/**
 * @brief Apply one post-step decision to integrator buffers.
 *
 * `positions` and `velocities` are the state after the step. `previous_*`
 * are the state before the step and are not modified. A NULL `fn` returns
 * `LAMMPC_BOUNDARY_ACCEPT` without reading the buffers. Accept restores
 * the post-step buffers if the callback wrote them. Reject copies
 * `previous_*` back with `memcpy`. Reflect leaves the callback's writes.
 * Any other code restores `previous_*` and returns -1.
 *
 * @return the decision, or -1 when the request or the callback code is
 * invalid. A NULL callback is valid and returns accept.
 */
int lammpc_hard_boundary_apply(LammpcHardBoundaryFn fn, void *user,
                               int n_atoms, double *positions,
                               double *velocities,
                               const double *previous_positions,
                               const double *previous_velocities,
                               const LammpcBoundaryBox *box);

#ifdef __cplusplus
}
#endif

#endif /* LAMMPC_BOUNDARY_H */
