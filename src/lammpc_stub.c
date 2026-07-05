/**
 * Standalone stub: the full lammpc ABI surface without a linked LAMMPS.
 * Every evaluation fails loudly; discovery (features, version, abi_version)
 * answers exactly like the embed build with available = 0.
 */
#include "lammpc.h"

#include <stdio.h>

#ifndef LAMMPC_VERSION_STRING
#define LAMMPC_VERSION_STRING "unknown"
#endif

static const char *g_last_error = "lammpc stub: LAMMPS library not linked";

const char *lammpc_last_error(void) { return g_last_error; }

static LammpcResult stub_fail(void) {
  LammpcResult r;
  r.ok = 0;
  r.energy_h = 0.0;
  snprintf(r.message, sizeof(r.message), "%s", g_last_error);
  return r;
}

int lammpc_set_params(const void *params_capnp,
                      size_t params_capnp_size_bytes) {
  (void)params_capnp;
  (void)params_capnp_size_bytes;
  return -1;
}

int lammpc_configure(const void *config_capnp,
                     size_t config_capnp_size_bytes) {
  (void)config_capnp;
  (void)config_capnp_size_bytes;
  return -1;
}

LammpcResult lammpc_energy_gradient(int n_atoms, const double *positions_ang,
                                    const int *atomic_numbers,
                                    const void *cell_ang, size_t cell_len,
                                    double *gradient_h_bohr) {
  (void)n_atoms;
  (void)positions_ang;
  (void)atomic_numbers;
  (void)cell_ang;
  (void)cell_len;
  (void)gradient_h_bohr;
  return stub_fail();
}

LammpcResult lammpc_energy(int n_atoms, const double *positions_ang,
                           const int *atomic_numbers, const void *cell_ang,
                           size_t cell_len) {
  (void)n_atoms;
  (void)positions_ang;
  (void)atomic_numbers;
  (void)cell_ang;
  (void)cell_len;
  return stub_fail();
}

LammpcResult lammpc_energy_forces(int n_atoms, const double *positions_ang,
                                  const int *atomic_numbers,
                                  const void *cell_ang, size_t cell_len,
                                  double *forces_h_bohr) {
  (void)n_atoms;
  (void)positions_ang;
  (void)atomic_numbers;
  (void)cell_ang;
  (void)cell_len;
  (void)forces_h_bohr;
  return stub_fail();
}

LammpcSession *lammpc_session_create(const void *params_capnp,
                                     size_t params_capnp_size_bytes) {
  (void)params_capnp;
  (void)params_capnp_size_bytes;
  return NULL;
}

int lammpc_session_set_params(LammpcSession *session, const void *params_capnp,
                              size_t params_capnp_size_bytes) {
  (void)session;
  (void)params_capnp;
  (void)params_capnp_size_bytes;
  return -1;
}

LammpcSession *lammpc_session_create_from_config(
    const void *config_capnp, size_t config_capnp_size_bytes) {
  (void)config_capnp;
  (void)config_capnp_size_bytes;
  return NULL;
}

int lammpc_session_configure(LammpcSession *session, const void *config_capnp,
                             size_t config_capnp_size_bytes) {
  (void)session;
  (void)config_capnp;
  (void)config_capnp_size_bytes;
  return -1;
}

void lammpc_session_destroy(LammpcSession *session) { (void)session; }

LammpcResult lammpc_session_calculate_result(
    LammpcSession *session, const void *force_input_capnp,
    size_t force_input_capnp_size_bytes, void *potential_result_capnp,
    size_t potential_result_capnp_capacity_bytes,
    size_t *potential_result_capnp_size_bytes) {
  (void)session;
  (void)force_input_capnp;
  (void)force_input_capnp_size_bytes;
  (void)potential_result_capnp;
  (void)potential_result_capnp_capacity_bytes;
  (void)potential_result_capnp_size_bytes;
  return stub_fail();
}

LammpcResult lammpc_calculate_result(const void *params_capnp,
                                     size_t params_capnp_size_bytes,
                                     const void *force_input_capnp,
                                     size_t force_input_capnp_size_bytes,
                                     void *potential_result_capnp,
                                     size_t potential_result_capnp_capacity_bytes,
                                     size_t *potential_result_capnp_size_bytes) {
  (void)params_capnp;
  (void)params_capnp_size_bytes;
  (void)force_input_capnp;
  (void)force_input_capnp_size_bytes;
  (void)potential_result_capnp;
  (void)potential_result_capnp_capacity_bytes;
  (void)potential_result_capnp_size_bytes;
  return stub_fail();
}

LammpcResult lammpc_calculate_result_from_config(
    const void *config_capnp, size_t config_capnp_size_bytes,
    const void *force_input_capnp, size_t force_input_capnp_size_bytes,
    void *potential_result_capnp,
    size_t potential_result_capnp_capacity_bytes,
    size_t *potential_result_capnp_size_bytes) {
  (void)config_capnp;
  (void)config_capnp_size_bytes;
  (void)force_input_capnp;
  (void)force_input_capnp_size_bytes;
  (void)potential_result_capnp;
  (void)potential_result_capnp_capacity_bytes;
  (void)potential_result_capnp_size_bytes;
  return stub_fail();
}

size_t lammpc_potential_result_size_for_force_input(
    const void *force_input_capnp, size_t force_input_capnp_size_bytes) {
  (void)force_input_capnp;
  (void)force_input_capnp_size_bytes;
  return 0;
}

const char *lammpc_version(void) { return "lammpc-stub/" LAMMPC_VERSION_STRING; }

int lammpc_abi_version(void) { return LAMMPC_ABI_VERSION; }

int lammpc_available(void) { return 0; }

void lammpc_finalize(void) {}
