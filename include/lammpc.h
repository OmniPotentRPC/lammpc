#ifndef LAMMPC_H
#define LAMMPC_H

#include "lammpc_features.h"

#include <stddef.h>

/**
 * @file lammpc.h
 * @brief Stable C ABI over LAMMPS as a dynlib (the lammpc shim).
 *
 * Matches the shared-library soversion and lammpc_abi_version(); bumps only on
 * incompatible signature changes. Conforms to the minimum potential ABI
 * profile (potentials-schema PROFILE.md): flat binary Cap'n Proto messages in
 * both directions, plain C calls in-process. The embed build calls the actual
 * liblammps functions directly -- instance creation, pair setup and pair
 * compute, position/force array access -- with no command-string parsing, no
 * input scripts, and no file traffic anywhere in the force loop.
 */
#define LAMMPC_ABI_VERSION 0

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Result of one evaluation; layout fixed by the ABI profile. */
typedef struct LammpcResult {
  int ok;            /**< 1 success, 0 failure */
  double energy_h;   /**< total energy in Hartree */
  char message[512]; /**< status / error, NUL-terminated */
} LammpcResult;

/** @brief Opaque persistent evaluation session (one LAMMPS instance). */
typedef struct LammpcSession LammpcSession;

/**
 * @brief Apply sticky interaction-model options from a flat `LammpsParams`
 * message (units/pair_style/pair_coeff/type map).
 *
 * @return 0 on success, non-zero if the embed shell is unavailable or the
 * message is invalid; details via `lammpc_last_error()`.
 */
int lammpc_set_params(const void *params_capnp,
                      size_t params_capnp_size_bytes);

/**
 * @brief Apply a full `PotentialConfig`: the `lammps` union arm wins where it
 * differs from schema defaults. A set `common` overlay (CommonMethodSpec)
 * has no classical-force-field lowering; any set overlay field fails loudly
 * via `lammpc_last_error()`.
 */
int lammpc_configure(const void *config_capnp,
                     size_t config_capnp_size_bytes);

LammpcResult lammpc_energy_gradient(int n_atoms, const double *positions_ang,
                                    const int *atomic_numbers,
                                    const void *cell_ang, size_t cell_len,
                                    double *gradient_h_bohr);

LammpcResult lammpc_energy(int n_atoms, const double *positions_ang,
                           const int *atomic_numbers, const void *cell_ang,
                           size_t cell_len);

LammpcResult lammpc_energy_forces(int n_atoms, const double *positions_ang,
                                  const int *atomic_numbers,
                                  const void *cell_ang, size_t cell_len,
                                  double *forces_h_bohr);

/** @brief Create a persistent session from flat `LammpsParams` bytes. */
LammpcSession *lammpc_session_create(const void *params_capnp,
                                     size_t params_capnp_size_bytes);

int lammpc_session_set_params(LammpcSession *session, const void *params_capnp,
                              size_t params_capnp_size_bytes);

/** @brief Create a session from a full `PotentialConfig` (see
 * `lammpc_configure()`). */
LammpcSession *lammpc_session_create_from_config(
    const void *config_capnp, size_t config_capnp_size_bytes);

int lammpc_session_configure(LammpcSession *session, const void *config_capnp,
                             size_t config_capnp_size_bytes);

void lammpc_session_destroy(LammpcSession *session);

/**
 * @brief One step on a persistent session: `ForceInput` bytes in,
 * `PotentialResult` bytes out. This is the BOMD force-loop entry point for
 * rgpot callers driving dynamics.
 */
LammpcResult lammpc_session_calculate_result(
    LammpcSession *session, const void *force_input_capnp,
    size_t force_input_capnp_size_bytes, void *potential_result_capnp,
    size_t potential_result_capnp_capacity_bytes,
    size_t *potential_result_capnp_size_bytes);

/** @brief One-shot params + step. */
LammpcResult lammpc_calculate_result(const void *params_capnp,
                                     size_t params_capnp_size_bytes,
                                     const void *force_input_capnp,
                                     size_t force_input_capnp_size_bytes,
                                     void *potential_result_capnp,
                                     size_t potential_result_capnp_capacity_bytes,
                                     size_t *potential_result_capnp_size_bytes);

/** @brief One-shot config + step (full `PotentialConfig` resolution). */
LammpcResult lammpc_calculate_result_from_config(
    const void *config_capnp, size_t config_capnp_size_bytes,
    const void *force_input_capnp, size_t force_input_capnp_size_bytes,
    void *potential_result_capnp,
    size_t potential_result_capnp_capacity_bytes,
    size_t *potential_result_capnp_size_bytes);

/**
 * @brief Byte count needed for a `PotentialResult` for the given `ForceInput`.
 *
 * Parses geometry only; does not create a LAMMPS instance. Returns 0 when the
 * message is invalid or too large for the C ABI.
 */
size_t lammpc_potential_result_size_for_force_input(
    const void *force_input_capnp, size_t force_input_capnp_size_bytes);

/**
 * @brief Write a Cap'n Proto `Capabilities` message describing this backend.
 *
 * Loaders negotiate against the message before dispatch: backend name and
 * version, ABI generation, availability, the calculate operations the ABI
 * serves, the `CommonMethodSpec` fields the overlay lowers (none for a
 * classical engine), and the `PotentialConfig` arms accepted. A stub build
 * reports the same operation surface with `available = false`.
 *
 * Returns 0 on success. On a too-small buffer (including the pure size query
 * `capabilities_capnp == NULL`, `capabilities_capnp_capacity_bytes == 0`)
 * returns -1 with `*capabilities_capnp_size_bytes` set to the required size.
 */
int lammpc_capabilities_result(void *capabilities_capnp,
                               size_t capabilities_capnp_capacity_bytes,
                               size_t *capabilities_capnp_size_bytes);

/** @brief Compiled library version string. */
const char *lammpc_version(void);

/**
 * @brief Thread-local diagnostic for the most recent int-returning
 * configuration call; empty string on success.
 */
const char *lammpc_last_error(void);

/** @brief Numeric ABI generation; compare against LAMMPC_ABI_VERSION. */
int lammpc_abi_version(void);

/** @brief 1 when the LAMMPS library is linked and usable. */
int lammpc_available(void);

/** @brief Release the owned LAMMPS runtime. */
void lammpc_finalize(void);

#ifdef __cplusplus
}
#endif

#endif /* LAMMPC_H */
