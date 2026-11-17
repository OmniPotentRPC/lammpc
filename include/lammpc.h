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
 * both directions, plain C calls in-process. The embed build uses direct
 * liblammps calls for atom arrays and pair evaluation. Setup still uses a
 * small command set for region/box/mass lines, formatted at `%.17g` for
 * floating values. Atom creation is a full-precision batched
 * `lammps_create_atoms` path with a hard atom-count consistency check.
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

/** @brief Coordinate convention used by a constrained MD plane. */
typedef enum LammpcPlaneCoordinates {
  /** @brief q_i = sqrt(m_i) x_i; the supplied normal is Euclidean in q. */
  LAMMPC_PLANE_MASS_WEIGHTED = 0,
  /** @brief Cartesian point and normal with the mass metric applied by MD. */
  LAMMPC_PLANE_CARTESIAN = 1
} LammpcPlaneCoordinates;

/**
 * @brief Controls for lammpc's generic projected BAOAB stepper.
 *
 * Metal units are used: time_step_ps in ps, thermal_energy_ev in eV, and
 * friction_per_ps in ps^-1. A zero seed is valid and deterministic.
 */
typedef struct LammpcMdConfig {
  double time_step_ps;
  double thermal_energy_ev;
  double friction_per_ps;
  unsigned long long seed;
  LammpcPlaneCoordinates coordinates;
} LammpcMdConfig;

/** @brief Scalar accounting for the current constrained MD state. */
typedef struct LammpcMdObservation {
  /** Force along the positive unit normal in the selected coordinates. */
  double constraint_force;
  /** Number of force evaluations since md_initialize. */
  size_t force_evaluations;
  /** Number of completed BAOAB steps. */
  size_t steps;
} LammpcMdObservation;

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
 * @brief Initialize constrained Langevin dynamics on a configured session.
 *
 * force_input supplies the atoms, positions, optional cell, and requested
 * eV/angstrom units. point and normal contain one value per Cartesian
 * coordinate. The initial position is restored to the plane, velocities are
 * drawn from the configured Maxwell distribution, and both velocity and
 * force are projected before returning.
 */
int lammpc_session_md_initialize(
    LammpcSession *session, const void *force_input_capnp,
    size_t force_input_capnp_size_bytes, const double *point,
    const double *normal, size_t coordinate_count,
    const LammpcMdConfig *config);

/** @brief Advance exactly step_count projected BAOAB steps. */
int lammpc_session_md_step(LammpcSession *session, size_t step_count);

/**
 * @brief Copy the current constrained MD state in Cartesian metal units.
 *
 * Every array has coordinate_count entries. Any output array may be NULL.
 * projected_forces_ev_ang contains the force used by the integrator.
 * removed_velocities_ang_ps contains the normal velocity removed by the most
 * recent projection.
 */
int lammpc_session_md_state(
    LammpcSession *session, double *positions_ang,
    double *velocities_ang_ps, double *projected_forces_ev_ang,
    double *removed_velocities_ang_ps, size_t coordinate_count,
    LammpcMdObservation *observation);

/**
 * @brief Evaluate raw arrays on a persistent session.
 *
 * Positions and the optional 3x3 row-major cell are in Angstrom. Energy is
 * returned in `LammpcResult::energy_h`, gradients in Hartree/Bohr, and the
 * optional row-major Cauchy stress in Hartree/Bohr^3. The atomic numbers
 * select the configured type map when the session first builds its LAMMPS
 * instance.
 */
LammpcResult lammpc_session_energy_gradient(
    LammpcSession *session, int n_atoms, const double *positions_ang,
    const int *atomic_numbers, const void *cell_ang, size_t cell_len,
    double *gradient_h_bohr, double *stress_h_bohr3);

/**
 * @brief One step on a persistent session: `ForceInput` bytes in,
 * `PotentialResult` bytes out. This is the BOMD force-loop entry point for
 * rgpot callers driving dynamics. Geometry setup fails loudly when LAMMPS
 * creates fewer or more atoms than requested.
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
