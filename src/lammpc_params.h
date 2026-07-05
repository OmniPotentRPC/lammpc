#pragma once

#include "Potentials.capnp.h"

#include <stddef.h>

#define LAMMPC_BOHR_TO_ANGSTROM 0.529177210903
#define LAMMPC_HARTREE_TO_EV 27.211386245988

/** Root a flat `ForceInput` message; caller releases via
 * `lammpc_params_release()`. */
int lammpc_force_input_root(const void *force_input_capnp,
                            size_t force_input_capnp_size_bytes,
                            struct capn *arena, ForceInput_ptr *force_input);

/** Root a flat `LammpsParams` message. */
int lammpc_params_root(const void *params_capnp,
                       size_t params_capnp_size_bytes, struct capn *arena,
                       LammpsParams_ptr *params);

void lammpc_params_release(struct capn *arena);

int lammpc_force_input_atom_count(ForceInput_ptr force_input, size_t *n_atoms,
                                  int *has_cell);

/** Copy geometry into caller buffers, converting lengths to Angstrom. */
int lammpc_force_input_copy_geometry(ForceInput_ptr force_input,
                                     double *positions_ang,
                                     int *atomic_numbers, size_t atom_capacity,
                                     double *cell_ang, int *has_cell);

/** Hartree->energyUnit and Hartree/Bohr->energyUnit/lengthUnit factors. */
int lammpc_force_input_result_factors(ForceInput_ptr force_input,
                                      double *energy_factor,
                                      double *force_factor);

size_t lammpc_potential_result_flat_size(size_t force_count);

/** Write energy + forces into a flat `PotentialResult` message. */
int lammpc_potential_result_write(double energy, const double *forces,
                                  size_t force_count,
                                  void *potential_result_capnp,
                                  size_t potential_result_capacity_bytes,
                                  size_t *potential_result_size_bytes);
