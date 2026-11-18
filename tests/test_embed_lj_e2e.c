/**
 * Real-liblammps E2E: an LJ dimer evaluated through the public ABI must
 * reproduce the analytic lj/cut energy and forces, including on a session
 * reused across geometry steps (the BOMD force-loop path).
 */
#include "lammpc.h"
#include "lammpc_params.h"

#include <math.h>
#include <setjmp.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include <cmocka.h>

static const double kEps = 0.2;   /* eV */
static const double kSigma = 2.5; /* Angstrom */

static double lj_energy(double r) {
  double sr6 = pow(kSigma / r, 6.0);
  return 4.0 * kEps * (sr6 * sr6 - sr6);
}

static double lj_force_magnitude(double r) {
  double sr6 = pow(kSigma / r, 6.0);
  return 24.0 * kEps * (2.0 * sr6 * sr6 - sr6) / r;
}

typedef struct MdBoundaryCounts {
  int calls;
  int reflections;
  int valid;
} MdBoundaryCounts;

typedef struct MdRejectBoundary {
  int calls;
  double proposals[5][6];
} MdRejectBoundary;

typedef struct MdShiftBoundary {
  int calls;
  double shift;
  double velocity_delta[6];
  double reflected_velocities[6];
} MdShiftBoundary;

static int reflect_md_boundary(void *user, int n_atoms, double *positions,
                               double *velocities,
                               const LammpcBoundaryBox *box) {
  MdBoundaryCounts *counts = (MdBoundaryCounts *)user;
  (void)positions;
  (void)velocities;
  counts->calls += 1;
  if (n_atoms != 2 || !(box->hi[0] > box->lo[0]))
    counts->valid = 0;
  counts->reflections += 1;
  return LAMMPC_BOUNDARY_REFLECT;
}

static int reject_md_boundary(void *user, int n_atoms, double *positions,
                              double *velocities,
                              const LammpcBoundaryBox *box) {
  MdRejectBoundary *boundary = (MdRejectBoundary *)user;
  (void)box;
  assert_int_equal(n_atoms, 2);
  if (boundary->calls < 5)
    memcpy(boundary->proposals[boundary->calls], positions,
           sizeof(boundary->proposals[0]));
  boundary->calls += 1;
  positions[0] += 100.0;
  velocities[0] += 100.0;
  return LAMMPC_BOUNDARY_REJECT;
}

static int invalid_md_boundary(void *user, int n_atoms, double *positions,
                               double *velocities,
                               const LammpcBoundaryBox *box) {
  MdRejectBoundary *boundary = (MdRejectBoundary *)user;
  (void)box;
  boundary->calls += 1;
  assert_int_equal(n_atoms, 2);
  positions[0] += 100.0;
  velocities[0] += 100.0;
  return 99;
}

static int shift_md_boundary(void *user, int n_atoms, double *positions,
                             double *velocities,
                             const LammpcBoundaryBox *box) {
  MdShiftBoundary *boundary = (MdShiftBoundary *)user;
  (void)box;
  boundary->calls += 1;
  assert_int_equal(n_atoms, 2);
  positions[0] -= boundary->shift;
  positions[3] += boundary->shift;
  for (int i = 0; i < 6; ++i) {
    velocities[i] += boundary->velocity_delta[i];
    boundary->reflected_velocities[i] = velocities[i];
  }
  return LAMMPC_BOUNDARY_REFLECT;
}

static size_t write_flat(struct capn *arena, capn_ptr root_struct,
                         unsigned char *dst, size_t dst_size) {
  capn_ptr root = capn_root(arena);
  assert_int_equal(capn_setp(root, 0, root_struct), 0);
  int written = capn_write_mem(arena, dst, dst_size, 0);
  assert_true(written > 0);
  return (size_t)written;
}

static size_t make_config(unsigned char *dst, size_t dst_size) {
  struct capn arena;
  capn_init_malloc(&arena);
  struct capn_segment *seg = capn_root(&arena).seg;

  LammpsParams_ptr params = new_LammpsParams(seg);
  capn_ptr coeffs = capn_new_ptr_list(seg, 1);
  capn_text coeff = {13, "1 1 0.2 2.5", 0};
  coeff.len = (int)strlen("1 1 0.2 2.5");
  assert_int_equal(capn_set_text(coeffs, 0, coeff), 0);
  capn_list32 type_map = capn_new_list32(seg, 1);
  capn_set32(type_map, 0, 29);

  struct LammpsParams pview;
  memset(&pview, 0, sizeof(pview));
  pview.unitsStyle = (capn_text){5, "metal", NULL};
  pview.atomStyle = (capn_text){6, "atomic", NULL};
  pview.pairStyle = (capn_text){11, "lj/cut 10.0", NULL};
  pview.pairCoeffs = coeffs;
  pview.typeToAtomicNumber = type_map;
  pview.newtonPair = 1;
  pview.boundary = (capn_text){5, "f f f", NULL};
  write_LammpsParams(&pview, params);

  PotentialConfig_ptr config = new_PotentialConfig(seg);
  struct PotentialConfig cview;
  memset(&cview, 0, sizeof(cview));
  cview.which = PotentialConfig_lammps;
  cview.lammps = params;
  write_PotentialConfig(&cview, config);

  size_t written = write_flat(&arena, config.p, dst, dst_size);
  capn_free(&arena);
  return written;
}

static size_t make_step_positions(const double coordinates[6],
                                  unsigned char *dst, size_t dst_size) {
  struct capn arena;
  capn_init_malloc(&arena);
  struct capn_segment *seg = capn_root(&arena).seg;

  ForceInput_ptr input = new_ForceInput(seg);
  capn_list64 pos = capn_new_list64(seg, 6);
  for (int i = 0; i < 6; ++i)
    capn_set64(pos, i, capn_from_f64(coordinates[i]));
  capn_list32 atmnrs = capn_new_list32(seg, 2);
  capn_set32(atmnrs, 0, 29);
  capn_set32(atmnrs, 1, 29);

  struct ForceInput view;
  memset(&view, 0, sizeof(view));
  view.pos = pos;
  view.atmnrs = atmnrs;
  view.lengthUnit = (capn_text){8, "angstrom", NULL};
  view.energyUnit = (capn_text){2, "eV", NULL};
  write_ForceInput(&view, input);

  size_t written = write_flat(&arena, input.p, dst, dst_size);
  capn_free(&arena);
  return written;
}

static size_t make_step(double separation_ang, unsigned char *dst,
                        size_t dst_size) {
  const double coordinates[6] = {0.0, 0.0, 0.0,
                                 separation_ang, 0.0, 0.0};
  return make_step_positions(coordinates, dst, dst_size);
}

static void decode_result(const unsigned char *buf, size_t size,
                          double *energy, double forces[6]) {
  struct capn arena;
  assert_int_equal(capn_init_mem(&arena, buf, size, 0), 0);
  PotentialResult_ptr result;
  result.p = capn_getp(capn_root(&arena), 0, 1);
  assert_int_equal(result.p.type, CAPN_STRUCT);
  struct PotentialResult view;
  read_PotentialResult(&view, result);
  *energy = view.energy;
  capn_resolve(&view.forces.p);
  assert_int_equal(view.forces.p.len, 6);
  for (int i = 0; i < 6; ++i)
    forces[i] = capn_to_f64(capn_get64(view.forces, i));
  capn_free(&arena);
}

static void check_step(LammpcSession *session, const unsigned char *config,
                       size_t config_size, double separation) {
  (void)config;
  (void)config_size;
  unsigned char step[4096];
  size_t step_size = make_step(separation, step, sizeof(step));

  size_t need = lammpc_potential_result_size_for_force_input(step, step_size);
  assert_true(need > 0);
  unsigned char *out = (unsigned char *)malloc(need);
  assert_non_null(out);
  size_t wrote = 0;
  LammpcResult r = lammpc_session_calculate_result(session, step, step_size,
                                                   out, need, &wrote);
  if (!r.ok)
    fail_msg("session step failed: %s", r.message);
  assert_true(wrote > 0);

  double energy = 0.0;
  double forces[6];
  decode_result(out, wrote, &energy, forces);
  free(out);

  double expected_e = lj_energy(separation);
  double expected_f = lj_force_magnitude(separation);
  assert_true(fabs(energy - expected_e) < 1e-9);
  /* Attractive at r > 2^(1/6) sigma: atom 2 pulled toward atom 1. */
  assert_true(fabs(forces[3] - expected_f) < 1e-9);
  assert_true(fabs(forces[0] + expected_f) < 1e-9);
  assert_true(fabs(forces[1]) < 1e-12 && fabs(forces[2]) < 1e-12);
  assert_true(fabs(forces[4]) < 1e-12 && fabs(forces[5]) < 1e-12);
}

static void test_lj_dimer_session_matches_analytic(void **state) {
  (void)state;
  assert_int_equal(lammpc_available(), 1);

  unsigned char config[8192];
  size_t config_size = make_config(config, sizeof(config));

  LammpcSession *session =
      lammpc_session_create_from_config(config, config_size);
  assert_non_null(session);

  check_step(session, config, config_size, 3.0);
  /* Same instance, new positions: the live per-step path. */
  check_step(session, config, config_size, 3.2);
  check_step(session, config, config_size, 2.9);

  lammpc_session_destroy(session);
}

static void test_lj_dimer_one_shot_config(void **state) {
  (void)state;
  unsigned char config[8192];
  size_t config_size = make_config(config, sizeof(config));
  unsigned char step[4096];
  size_t step_size = make_step(3.0, step, sizeof(step));

  size_t need = lammpc_potential_result_size_for_force_input(step, step_size);
  assert_true(need > 0);
  unsigned char *out = (unsigned char *)malloc(need);
  assert_non_null(out);
  size_t wrote = 0;
  LammpcResult r = lammpc_calculate_result_from_config(
      config, config_size, step, step_size, out, need, &wrote);
  if (!r.ok)
    fail_msg("one-shot failed: %s", r.message);

  double energy = 0.0;
  double forces[6];
  decode_result(out, wrote, &energy, forces);
  free(out);
  assert_true(fabs(energy - lj_energy(3.0)) < 1e-9);
}

static void test_projected_md_stays_on_plane(void **state) {
  (void)state;
  unsigned char config[8192];
  size_t config_size = make_config(config, sizeof(config));
  unsigned char step[4096];
  size_t step_size = make_step(3.0, step, sizeof(step));
  LammpcSession *session =
      lammpc_session_create_from_config(config, config_size);
  assert_non_null(session);

  const double inverse_sqrt_two = 1.0 / sqrt(2.0);
  const double point[6] = {0.0, 0.0, 0.0, 3.0, 0.0, 0.0};
  const double normal[6] = {-inverse_sqrt_two, 0.0, 0.0,
                            inverse_sqrt_two,  0.0, 0.0};
  const LammpcMdConfig controls = {
      .time_step_ps = 0.0002,
      .thermal_energy_ev = 0.01,
      .friction_per_ps = 10.0,
      .seed = 17,
      .coordinates = LAMMPC_PLANE_MASS_WEIGHTED,
  };
  assert_int_equal(lammpc_session_md_initialize(
                       session, step, step_size, point, normal, 6, &controls),
                   0);
  assert_int_equal(lammpc_session_md_step(session, 200), 0);

  double positions[6];
  double forces[6];
  LammpcMdObservation observation;
  assert_int_equal(lammpc_session_md_state(session, positions, NULL, forces,
                                           NULL, 6, &observation),
                   0);
  double distance = 0.0;
  double projected_force = 0.0;
  for (int i = 0; i < 6; ++i) {
    distance += normal[i] * (positions[i] - point[i]);
    projected_force += normal[i] * forces[i];
  }
  assert_true(fabs(distance) < 1e-12);
  assert_true(fabs(projected_force) < 1e-12);
  assert_int_equal(observation.steps, 200);
  assert_int_equal(observation.force_evaluations, 201);
  lammpc_session_destroy(session);
}

static void test_md_step_applies_hard_boundary_each_step(void **state) {
  (void)state;
  unsigned char config[8192];
  size_t config_size = make_config(config, sizeof(config));
  unsigned char step[4096];
  size_t step_size = make_step(3.0, step, sizeof(step));
  LammpcSession *session =
      lammpc_session_create_from_config(config, config_size);
  assert_non_null(session);

  const double inverse_sqrt_two = 1.0 / sqrt(2.0);
  const double point[6] = {0.0, 0.0, 0.0, 3.0, 0.0, 0.0};
  const double normal[6] = {-inverse_sqrt_two, 0.0, 0.0,
                            inverse_sqrt_two,  0.0, 0.0};
  const LammpcMdConfig controls = {
      .time_step_ps = 0.0002,
      .thermal_energy_ev = 0.01,
      .friction_per_ps = 10.0,
      .seed = 19,
      .coordinates = LAMMPC_PLANE_MASS_WEIGHTED,
  };
  MdBoundaryCounts counts = {.calls = 0, .reflections = 0, .valid = 1};
  const int steps = 7;

  assert_int_equal(lammpc_session_md_initialize(
                       session, step, step_size, point, normal, 6, &controls),
                   0);
  assert_int_equal(lammpc_session_set_hard_boundary(
                       session, reflect_md_boundary, &counts),
                   0);
  for (int step_index = 0; step_index < steps; ++step_index)
    assert_int_equal(lammpc_session_md_step(session, 1), 0);

  LammpcMdObservation observation;
  assert_int_equal(
      lammpc_session_md_state(session, NULL, NULL, NULL, NULL, 6, &observation),
      0);
  assert_int_equal(counts.calls, steps);
  assert_int_equal(counts.reflections, steps);
  assert_int_equal(counts.valid, 1);
  assert_int_equal(observation.steps, (size_t)steps);
  assert_int_equal(observation.force_evaluations, 2u * (size_t)steps + 1u);
  lammpc_session_destroy(session);
}

static void test_md_reject_restores_pre_step_force_state(void **state) {
  (void)state;
  unsigned char config[8192];
  unsigned char step[4096];
  const size_t config_size = make_config(config, sizeof(config));
  const size_t step_size = make_step(3.0, step, sizeof(step));
  const double inverse_sqrt_two = 1.0 / sqrt(2.0);
  const double point[6] = {0.0, 0.0, 0.0, 3.0, 0.0, 0.0};
  const double normal[6] = {-inverse_sqrt_two, 0.0, 0.0,
                            inverse_sqrt_two,  0.0, 0.0};
  const LammpcMdConfig controls = {
      .time_step_ps = 0.0002,
      .thermal_energy_ev = 0.01,
      .friction_per_ps = 10.0,
      .seed = 23,
      .coordinates = LAMMPC_PLANE_MASS_WEIGHTED,
  };
  LammpcSession *rejected =
      lammpc_session_create_from_config(config, config_size);
  LammpcSession *reference =
      lammpc_session_create_from_config(config, config_size);
  assert_non_null(rejected);
  assert_non_null(reference);
  assert_int_equal(lammpc_session_md_initialize(rejected, step, step_size,
                                                 point, normal, 6, &controls),
                   0);
  assert_int_equal(lammpc_session_md_initialize(reference, step, step_size,
                                                 point, normal, 6, &controls),
                   0);

  MdRejectBoundary boundary = {.calls = 0, .proposals = {{0.0}}};
  assert_int_equal(lammpc_session_set_hard_boundary(
                       rejected, reject_md_boundary, &boundary),
                   0);
  assert_int_equal(lammpc_session_md_step(rejected, 5), 0);

  double rejected_positions[6], rejected_velocities[6], rejected_forces[6];
  double rejected_removed[6], reference_positions[6], reference_velocities[6];
  double reference_forces[6], reference_removed[6];
  LammpcMdObservation rejected_observation, reference_observation;
  assert_int_equal(lammpc_session_md_state(
                       rejected, rejected_positions, rejected_velocities,
                       rejected_forces, rejected_removed, 6,
                       &rejected_observation),
                   0);
  assert_int_equal(lammpc_session_md_state(
                       reference, reference_positions, reference_velocities,
                       reference_forces, reference_removed, 6,
                       &reference_observation),
                   0);
  assert_int_equal(boundary.calls, 5);
  for (int proposal = 1; proposal < boundary.calls; ++proposal)
    assert_memory_not_equal(boundary.proposals[proposal - 1],
                            boundary.proposals[proposal],
                            sizeof(boundary.proposals[proposal]));
  for (int i = 0; i < 6; ++i) {
    assert_true(fabs(rejected_positions[i] - reference_positions[i]) < 1e-14);
    assert_true(fabs(rejected_velocities[i] - reference_velocities[i]) < 1e-14);
    assert_true(fabs(rejected_forces[i] - reference_forces[i]) < 1e-14);
    assert_true(fabs(rejected_removed[i] - reference_removed[i]) < 1e-14);
  }
  assert_true(fabs(rejected_observation.constraint_force -
                   reference_observation.constraint_force) < 1e-14);
  assert_int_equal(rejected_observation.force_evaluations,
                   reference_observation.force_evaluations + 5u);
  assert_int_equal(rejected_observation.steps, 5u);
  lammpc_session_destroy(rejected);
  lammpc_session_destroy(reference);
}

static void test_md_invalid_decision_restores_pre_step_state(void **state) {
  (void)state;
  unsigned char config[8192];
  unsigned char step[4096];
  const size_t config_size = make_config(config, sizeof(config));
  const size_t step_size = make_step(3.0, step, sizeof(step));
  const double inverse_sqrt_two = 1.0 / sqrt(2.0);
  const double point[6] = {0.0, 0.0, 0.0, 3.0, 0.0, 0.0};
  const double normal[6] = {-inverse_sqrt_two, 0.0, 0.0,
                            inverse_sqrt_two,  0.0, 0.0};
  const LammpcMdConfig controls = {
      .time_step_ps = 0.0002,
      .thermal_energy_ev = 0.01,
      .friction_per_ps = 10.0,
      .seed = 31,
      .coordinates = LAMMPC_PLANE_MASS_WEIGHTED,
  };
  LammpcSession *invalid =
      lammpc_session_create_from_config(config, config_size);
  LammpcSession *reference =
      lammpc_session_create_from_config(config, config_size);
  assert_non_null(invalid);
  assert_non_null(reference);
  assert_int_equal(lammpc_session_md_initialize(invalid, step, step_size, point,
                                                 normal, 6, &controls),
                   0);
  assert_int_equal(lammpc_session_md_initialize(reference, step, step_size,
                                                 point, normal, 6, &controls),
                   0);

  MdRejectBoundary boundary = {.calls = 0, .proposals = {{0.0}}};
  assert_int_equal(lammpc_session_set_hard_boundary(
                       invalid, invalid_md_boundary, &boundary),
                   0);
  assert_int_equal(lammpc_session_md_step(invalid, 1), -1);
  assert_non_null(strstr(lammpc_last_error(), "hard boundary callback failed"));

  double invalid_positions[6], invalid_velocities[6], invalid_forces[6];
  double invalid_removed[6], reference_positions[6], reference_velocities[6];
  double reference_forces[6], reference_removed[6];
  LammpcMdObservation invalid_observation, reference_observation;
  assert_int_equal(lammpc_session_md_state(
                       invalid, invalid_positions, invalid_velocities,
                       invalid_forces, invalid_removed, 6,
                       &invalid_observation),
                   0);
  assert_int_equal(lammpc_session_md_state(
                       reference, reference_positions, reference_velocities,
                       reference_forces, reference_removed, 6,
                       &reference_observation),
                   0);
  assert_int_equal(boundary.calls, 1);
  for (int i = 0; i < 6; ++i) {
    assert_memory_equal(&invalid_positions[i], &reference_positions[i],
                        sizeof(double));
    assert_memory_equal(&invalid_velocities[i], &reference_velocities[i],
                        sizeof(double));
    assert_memory_equal(&invalid_forces[i], &reference_forces[i],
                        sizeof(double));
    assert_memory_equal(&invalid_removed[i], &reference_removed[i],
                        sizeof(double));
  }
  assert_memory_equal(&invalid_observation.constraint_force,
                      &reference_observation.constraint_force, sizeof(double));
  assert_int_equal(invalid_observation.force_evaluations,
                   reference_observation.force_evaluations + 1u);
  assert_int_equal(invalid_observation.steps, reference_observation.steps);
  lammpc_session_destroy(invalid);
  lammpc_session_destroy(reference);
}

static void test_md_reflect_reprojects_and_reevaluates(void **state) {
  (void)state;
  unsigned char config[8192];
  unsigned char step[4096];
  const size_t config_size = make_config(config, sizeof(config));
  const size_t step_size = make_step(3.0, step, sizeof(step));
  const double inverse_sqrt_two = 1.0 / sqrt(2.0);
  const double point[6] = {0.0, 0.0, 0.0, 3.0, 0.0, 0.0};
  const double normal[6] = {-inverse_sqrt_two, 0.0, 0.0,
                            inverse_sqrt_two,  0.0, 0.0};
  const LammpcMdConfig controls = {
      .time_step_ps = 0.0002,
      .thermal_energy_ev = 0.01,
      .friction_per_ps = 10.0,
      .seed = 29,
      .coordinates = LAMMPC_PLANE_MASS_WEIGHTED,
  };
  LammpcSession *reflected =
      lammpc_session_create_from_config(config, config_size);
  assert_non_null(reflected);
  assert_int_equal(lammpc_session_md_initialize(reflected, step, step_size,
                                                 point, normal, 6, &controls),
                   0);
  MdShiftBoundary boundary = {
      .calls = 0,
      .shift = 0.05,
      .velocity_delta = {2.0, -0.5, 0.25, -1.0, 0.75, -0.125},
      .reflected_velocities = {0.0},
  };
  assert_int_equal(lammpc_session_set_hard_boundary(
                       reflected, shift_md_boundary, &boundary),
                   0);
  assert_int_equal(lammpc_session_md_step(reflected, 1), 0);

  double positions[6], velocities[6], forces[6], removed[6];
  LammpcMdObservation observation;
  assert_int_equal(lammpc_session_md_state(reflected, positions, velocities,
                                           forces, removed, 6, &observation),
                   0);
  double distance = 0.0;
  double observed_kinetic_energy = 0.0;
  double expected_kinetic_energy = 0.0;
  double normal_component = 0.0;
  for (int i = 0; i < 6; ++i)
    normal_component += normal[i] * boundary.reflected_velocities[i];
  for (int i = 0; i < 6; ++i) {
    const double expected_velocity =
        boundary.reflected_velocities[i] - normal[i] * normal_component;
    const double expected_removed = normal[i] * normal_component;
    distance += normal[i] * (positions[i] - point[i]);
    observed_kinetic_energy += 0.5 * 63.546 * 1.0364269656262175e-4 *
                               velocities[i] * velocities[i];
    expected_kinetic_energy += 0.5 * 63.546 * 1.0364269656262175e-4 *
                               expected_velocity * expected_velocity;
    assert_true(fabs(velocities[i] - expected_velocity) < 1e-12);
    assert_true(fabs(removed[i] - expected_removed) < 1e-12);
  }
  assert_int_equal(boundary.calls, 1);
  assert_true(fabs(distance) < 1e-12);
  assert_true(fabs(observed_kinetic_energy - expected_kinetic_energy) < 1e-12);
  assert_int_equal(observation.force_evaluations, 3);

  unsigned char reflected_step[4096];
  const size_t reflected_step_size =
      make_step_positions(positions, reflected_step, sizeof(reflected_step));
  LammpcSession *fresh =
      lammpc_session_create_from_config(config, config_size);
  assert_non_null(fresh);
  assert_int_equal(lammpc_session_md_initialize(
                       fresh, reflected_step, reflected_step_size, positions,
                       normal, 6, &controls),
                   0);
  double fresh_forces[6];
  assert_int_equal(lammpc_session_md_state(fresh, NULL, NULL, fresh_forces,
                                           NULL, 6, NULL),
                   0);
  for (int i = 0; i < 6; ++i)
    assert_true(fabs(forces[i] - fresh_forces[i]) < 1e-12);
  lammpc_session_destroy(fresh);
  lammpc_session_destroy(reflected);
}

static void test_overlay_rejected(void **state) {
  (void)state;
  struct capn arena;
  capn_init_malloc(&arena);
  struct capn_segment *seg = capn_root(&arena).seg;

  LammpsParams_ptr params = new_LammpsParams(seg);
  struct LammpsParams pview;
  memset(&pview, 0, sizeof(pview));
  pview.pairStyle = (capn_text){11, "lj/cut 10.0", NULL};
  capn_list32 type_map = capn_new_list32(seg, 1);
  capn_set32(type_map, 0, 29);
  pview.typeToAtomicNumber = type_map;
  write_LammpsParams(&pview, params);

  CommonMethodSpec_ptr common = new_CommonMethodSpec(seg);
  struct CommonMethodSpec cm;
  memset(&cm, 0, sizeof(cm));
  cm.basisSet = (capn_text){9, "planewave", NULL};
  write_CommonMethodSpec(&cm, common);

  PotentialConfig_ptr config = new_PotentialConfig(seg);
  struct PotentialConfig cview;
  memset(&cview, 0, sizeof(cview));
  cview.which = PotentialConfig_lammps;
  cview.lammps = params;
  cview.common = common;
  write_PotentialConfig(&cview, config);

  unsigned char buf[8192];
  size_t written = write_flat(&arena, config.p, buf, sizeof(buf));
  capn_free(&arena);

  assert_null(lammpc_session_create_from_config(buf, written));
  assert_non_null(strstr(lammpc_last_error(), "overlay"));
  assert_int_not_equal(lammpc_configure(buf, written), 0);
}

int main(void) {
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_lj_dimer_session_matches_analytic),
      cmocka_unit_test(test_lj_dimer_one_shot_config),
      cmocka_unit_test(test_projected_md_stays_on_plane),
      cmocka_unit_test(test_md_step_applies_hard_boundary_each_step),
      cmocka_unit_test(test_md_reject_restores_pre_step_force_state),
      cmocka_unit_test(test_md_invalid_decision_restores_pre_step_state),
      cmocka_unit_test(test_md_reflect_reprojects_and_reevaluates),
      cmocka_unit_test(test_overlay_rejected),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
