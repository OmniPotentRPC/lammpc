#ifndef LAMMPC_H
#define LAMMPC_H

#include "lammpc_features.h"
#include "lammpc_boundary.h"

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
/** Numeric ABI generation used by this header and shared library. */
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
  double time_step_ps;       /**< Finite positive time step in ps. */
  double thermal_energy_ev;  /**< Finite positive thermal energy in eV. */
  double friction_per_ps;    /**< Finite positive friction in ps^-1. */
  unsigned long long seed;   /**< Deterministic RNG seed. Zero is valid. */
  LammpcPlaneCoordinates coordinates; /**< Plane coordinate convention. */
} LammpcMdConfig;

/** @brief Scalar accounting for the current constrained MD state. */
typedef struct LammpcMdObservation {
  /** Force along the positive unit normal in the selected coordinates. */
  double constraint_force;
  /** Force evaluations since md_initialize, including rejected/failed
   * attempts and reflected-state re-evaluations. */
  size_t force_evaluations;
  /** Successful md_step iterations; boundary-rejected iterations count. */
  size_t steps;
} LammpcMdObservation;

/**
 * @brief Unconstrained NVE observables in metal units.
 *
 * potential_energy_ev and kinetic_energy_ev are read from LAMMPS thermo
 * after the most recent NVE propagation.
 */
typedef struct LammpcNveObservation {
  double potential_energy_ev; /**< Potential energy in eV. */
  double kinetic_energy_ev;   /**< Kinetic energy in eV. */
} LammpcNveObservation;

/**
 * @brief Legacy global-configuration entry point.
 *
 * The current ABI uses session-scoped configuration only. This function keeps
 * its symbol for discovery. Both builds ignore both arguments and return -1.
 * The embed build stores a message saying that global parameters are not
 * supported in `lammpc_last_error()`. The stub build does not modify
 * `lammpc_last_error()`.
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

/**
 * @brief Legacy raw-array gradient entry point.
 *
 * Positions use angstroms and gradients use Hartree/Bohr. Both the embed and
 * stub builds return `ok = 0` (the embed requires a session and stores an error
 * in `LammpcResult.message` and `lammpc_last_error()`, while the stub copies
 * its static unavailable message into `LammpcResult.message`); use
 * `lammpc_session_calculate_result()` instead.
 */
LammpcResult lammpc_energy_gradient(int n_atoms, const double *positions_ang,
                                    const int *atomic_numbers,
                                    const void *cell_ang, size_t cell_len,
                                    double *gradient_h_bohr);

/**
 * @brief Legacy raw-array energy entry point.
 *
 * Both the embed and stub builds return `ok = 0`; use
 * `lammpc_session_calculate_result()` instead.
 */
LammpcResult lammpc_energy(int n_atoms, const double *positions_ang,
                           const int *atomic_numbers, const void *cell_ang,
                           size_t cell_len);

/**
 * @brief Legacy raw-array force entry point.
 *
 * Positions use angstroms and forces use Hartree/Bohr. Both the embed and stub
 * builds return `ok = 0`; use `lammpc_session_calculate_result()` instead.
 */
LammpcResult lammpc_energy_forces(int n_atoms, const double *positions_ang,
                                  const int *atomic_numbers,
                                  const void *cell_ang, size_t cell_len,
                                  double *forces_h_bohr);

/** @brief Create a persistent session from flat `LammpsParams` bytes. */
LammpcSession *lammpc_session_create(const void *params_capnp,
                                     size_t params_capnp_size_bytes);

/**
 * @brief Replace session parameters and discard any live geometry.
 *
 * @return 0 on success or -1 for a null session, invalid message, unsupported
 * units or atom style, or missing required fields.
 */
int lammpc_session_set_params(LammpcSession *session, const void *params_capnp,
                              size_t params_capnp_size_bytes);

/** @brief Create a session from a full `PotentialConfig` (see
 * `lammpc_configure()`). */
LammpcSession *lammpc_session_create_from_config(
    const void *config_capnp, size_t config_capnp_size_bytes);

/**
 * @brief Replace session configuration and discard any live geometry.
 *
 * @return 0 on success or -1 when the session or configuration is invalid.
 */
int lammpc_session_configure(LammpcSession *session, const void *config_capnp,
                             size_t config_capnp_size_bytes);

/** @brief Release a session and its LAMMPS instance. NULL is accepted. */
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

/**
 * @brief Advance exactly step_count projected BAOAB steps.
 *
 * When a hard-boundary callback is armed, it receives the unwrapped
 * post-step MD positions once per step. A reflected position is projected
 * back to the configured plane and its force is re-evaluated. A rejected
 * iteration restores the accepted state but retains the consumed Langevin
 * draws, and it counts as a successful md_step iteration.
 */
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
 * @brief Set unconstrained NVE positions for an existing session geometry.
 *
 * coordinate_count must match 3 * atom_count of the live session geometry.
 * If cell_ang is non-NULL then cell_len must be 9 and the simulation box is
 * updated before positions are scattered. When the session has a cell, far
 * periodic images are folded through the full lattice (h_inv) along periodic
 * axes before scattering, while non-periodic axes are left untouched.
 */
int lammpc_session_nve_set_positions(LammpcSession *session,
                                     const double *positions_ang,
                                     size_t coordinate_count,
                                     const double *cell_ang, size_t cell_len);

/**
 * @brief Set unconstrained NVE Cartesian velocities in angstrom/ps.
 *
 * coordinate_count must match 3 * atom_count of the live session geometry.
 */
int lammpc_session_nve_set_velocities(LammpcSession *session,
                                      const double *velocities_ang_ps,
                                      size_t coordinate_count);

/**
 * @brief Advance unconstrained NVE dynamics with LAMMPS's own integrator.
 *
 * Runs step_count velocity-Verlet steps (`fix nve`) at the given timestep
 * in ps. No projection or position restoration is applied. This entry point
 * ignores a hard-boundary callback armed on the session.
 */
int lammpc_session_nve_step(LammpcSession *session, size_t step_count,
                            double time_step_ps);

/**
 * @brief Copy unconstrained NVE state in Cartesian metal units.
 *
 * positions_ang and velocities_ang_ps each have coordinate_count entries.
 * Any output pointer may be NULL.
 */
int lammpc_session_nve_state(LammpcSession *session, double *positions_ang,
                             double *velocities_ang_ps,
                             size_t coordinate_count,
                             LammpcNveObservation *observation);

/**
 * @brief One step on a persistent session: `ForceInput` bytes in,
 * `PotentialResult` bytes out. This is the BOMD force-loop entry point for
 * rgpot callers driving dynamics. Geometry setup fails loudly when LAMMPS
 * creates fewer or more atoms than requested. An empty
 * `ForceInput.energyUnit` means electronvolt (`eV`).
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
 * In the embed build, parses geometry only without creating a LAMMPS instance
 * and returns 0 when the message is invalid or too large for the C ABI. In the
 * stub build, always returns 0. Does not set or clear `lammpc_last_error()`.
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
 * If allocating or encoding the message fails, it also returns -1 and leaves
 * `*capabilities_capnp_size_bytes` at 0. It returns -1 without writing anything
 * when `capabilities_capnp_size_bytes` is `NULL`.
 */
int lammpc_capabilities_result(void *capabilities_capnp,
                               size_t capabilities_capnp_capacity_bytes,
                               size_t *capabilities_capnp_size_bytes);

/** @brief Compiled library version string. */
const char *lammpc_version(void);

/**
 * @brief Diagnostic string for the calling thread's most recent error-updating
 * API call.
 *
 * In the embed build, returns a non-NULL, NUL-terminated pointer into a
 * 512-byte thread-local buffer (messages longer than 511 bytes are truncated).
 * The contents are overwritten or cleared by the next error-updating call on
 * the same thread, and the buffer lifetime ends when that thread exits; callers
 * must copy the string if needed and must not free or modify it.
 *
 * In the embed build, session creation (`lammpc_session_create()`,
 * `lammpc_session_create_from_config()`), configuration
 * (`lammpc_configure()`, `lammpc_session_set_params()`,
 * `lammpc_session_configure()`), session and one-shot evaluation
 * (`lammpc_session_calculate_result()`, `lammpc_calculate_result()`,
 * `lammpc_calculate_result_from_config()`), MD, NVE, hard-boundary hook, and
 * `lammpc_session_integrate()` calls store a message here on failure and clear
 * it to an empty string on success.
 *
 * `lammpc_set_params()` and the legacy raw-array calls (`lammpc_energy()`,
 * `lammpc_energy_forces()`, `lammpc_energy_gradient()`) always fail in both
 * builds. In the embed build, each stores an error message here and never
 * clears it; the raw-array calls also copy the message into
 * `LammpcResult.message`. In the stub build, `lammpc_set_params()` returns -1
 * and produces no `LammpcResult`, while `lammpc_energy()`,
 * `lammpc_energy_forces()` and `lammpc_energy_gradient()` return `ok = 0` with
 * `LammpcResult.message` copied from the static unavailable string. Neither
 * stub path sets or clears `lammpc_last_error()`.
 *
 * `lammpc_hard_boundary_apply()`, `lammpc_capabilities_result()`,
 * `lammpc_abi_version()`, `lammpc_available()`,
 * `lammpc_potential_result_size_for_force_input()`, and
 * `lammpc_session_destroy()` neither set nor clear this buffer. The same
 * applies to `lammpc_finalize()`, `lammpc_version()`, `lammpc_last_error()`,
 * and the feature table queries (`lammpc_feature_count()`,
 * `lammpc_feature_table()`, `lammpc_feature_find()`).
 *
 * In the stub build, `lammpc_last_error()` always returns a fixed read-only
 * static unavailable string, and no call sets or clears it.
 */
const char *lammpc_last_error(void);

/** @brief Numeric ABI generation; compare against LAMMPC_ABI_VERSION. */
int lammpc_abi_version(void);

/** @brief 1 when the LAMMPS library is linked and usable. */
int lammpc_available(void);

/** @brief Release the owned LAMMPS runtime. */
void lammpc_finalize(void);

/**
 * @brief Arm a host hard-boundary callback on one session.
 *
 * The callback stays dormant until `lammpc_session_integrate()` or
 * `lammpc_session_md_step()` advances the session. Passing a NULL function
 * clears it. Clearing removes the embed fix when the LAMMPS instance is
 * already live, so a later force evaluation does not carry it. The force
 * loop never calls the callback.
 *
 * @return 0 on success, -1 when `session` is NULL or the embed fix cannot
 * be removed.
 */
int lammpc_session_set_hard_boundary(LammpcSession *session,
                                     LammpcHardBoundaryFn fn, void *user);

/** @brief Clear a hard-boundary callback. Equivalent to setting NULL. */
int lammpc_session_clear_hard_boundary(LammpcSession *session);

/**
 * @brief Run `nsteps` LAMMPS integration steps on a live session.
 *
 * The instance must already exist and must already have a time integrator
 * and a timestep (typically `fix nve` and `timestep` in extra setup).
 * When no callback is armed this is only that run: no fix is added and no
 * coordinate snapshot is taken. When a callback is armed, fix
 * `lammpc/boundary` is installed for this run and removed before return.
 * It calls the host after every step.
 *
 * @return 0 on success, -1 on a missing instance, a negative step count,
 * or a LAMMPS error (see `lammpc_last_error()`).
 */
int lammpc_session_integrate(LammpcSession *session, int nsteps);

#ifdef __cplusplus
}
#endif

#endif /* LAMMPC_H */
