#include "lammpc.h"
#include "lammpc_params.h"

#include <math.h>
#include <setjmp.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cmath>
#include <vector>

/*
 * Reuse the production embed implementation directly so this test can inspect
 * the live LAMMPS handle and assert post-create coordinates bit-for-bit.
 */
#include "../src/lammpc_embed_api.cpp"

extern "C" {
#include <cmocka.h>
}

static size_t write_flat(struct capn *arena, capn_ptr root_struct,
                         unsigned char *dst, size_t dst_size) {
  capn_ptr root = capn_root(arena);
  assert_int_equal(capn_setp(root, 0, root_struct), 0);
  int written = capn_write_mem(arena, dst, dst_size, 0);
  assert_true(written > 0);
  return (size_t)written;
}

static size_t make_config(const char *boundary, unsigned char *dst,
                          size_t dst_size) {
  struct capn arena;
  capn_init_malloc(&arena);
  struct capn_segment *seg = capn_root(&arena).seg;

  LammpsParams_ptr params = new_LammpsParams(seg);
  capn_ptr coeffs = capn_new_ptr_list(seg, 1);
  capn_text coeff = {3, "* *", 0};
  assert_int_equal(capn_set_text(coeffs, 0, coeff), 0);
  capn_list32 type_map = capn_new_list32(seg, 1);
  capn_set32(type_map, 0, 14);
  capn_list64 masses = capn_new_list64(seg, 1);
  capn_set64(masses, 0, capn_from_f64(28.0855));

  struct LammpsParams pview;
  memset(&pview, 0, sizeof(pview));
  pview.unitsStyle = (capn_text){5, "metal", NULL};
  pview.atomStyle = (capn_text){6, "atomic", NULL};
  pview.pairStyle = (capn_text){8, "zero 8.0", NULL};
  pview.pairCoeffs = coeffs;
  pview.typeToAtomicNumber = type_map;
  pview.masses = masses;
  pview.newtonPair = 1;
  pview.boundary = (capn_text){(int)strlen(boundary), (char *)boundary, NULL};
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

static size_t make_lj_config(unsigned char *dst, size_t dst_size) {
  struct capn arena;
  capn_init_malloc(&arena);
  struct capn_segment *seg = capn_root(&arena).seg;

  LammpsParams_ptr params = new_LammpsParams(seg);
  capn_ptr coeffs = capn_new_ptr_list(seg, 1);
  const char *coeff_value = "1 1 0.2 2.5";
  capn_text coeff = {(int)strlen(coeff_value), (char *)coeff_value, 0};
  assert_int_equal(capn_set_text(coeffs, 0, coeff), 0);
  capn_list32 type_map = capn_new_list32(seg, 1);
  capn_set32(type_map, 0, 29);
  capn_list64 masses = capn_new_list64(seg, 1);
  capn_set64(masses, 0, capn_from_f64(63.546));

  struct LammpsParams pview;
  memset(&pview, 0, sizeof(pview));
  pview.unitsStyle = (capn_text){5, "metal", NULL};
  pview.atomStyle = (capn_text){6, "atomic", NULL};
  pview.pairStyle = (capn_text){11, "lj/cut 10.0", NULL};
  pview.pairCoeffs = coeffs;
  pview.typeToAtomicNumber = type_map;
  pview.masses = masses;
  pview.newtonPair = 1;
  pview.boundary = (capn_text){5, (char *)"f f f", NULL};
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

static size_t make_step(const double *positions_ang, const int *z, int n_atoms,
                        const double *box9, unsigned char *dst,
                        size_t dst_size) {
  struct capn arena;
  capn_init_malloc(&arena);
  struct capn_segment *seg = capn_root(&arena).seg;

  ForceInput_ptr input = new_ForceInput(seg);
  capn_list64 pos = capn_new_list64(seg, n_atoms * 3);
  for (int i = 0; i < n_atoms * 3; ++i)
    capn_set64(pos, i, capn_from_f64(positions_ang[i]));
  capn_list32 atmnrs = capn_new_list32(seg, n_atoms);
  for (int i = 0; i < n_atoms; ++i)
    capn_set32(atmnrs, i, z[i]);

  struct ForceInput view;
  memset(&view, 0, sizeof(view));
  view.pos = pos;
  view.atmnrs = atmnrs;
  capn_list64 box = capn_new_list64(seg, 9);
  for (int i = 0; i < 9; ++i)
    capn_set64(box, i, capn_from_f64(box9[i]));
  view.box = box;
  view.lengthUnit = (capn_text){8, "angstrom", NULL};
  view.energyUnit = (capn_text){2, "eV", NULL};
  write_ForceInput(&view, input);

  size_t written = write_flat(&arena, input.p, dst, dst_size);
  capn_free(&arena);
  return written;
}

static LammpcSession *session_from_config(const char *boundary) {
  unsigned char config[16384];
  size_t config_size = make_config(boundary, config, sizeof(config));
  LammpcSession *session =
      lammpc_session_create_from_config(config, config_size);
  assert_non_null(session);
  return session;
}

static LammpcResult run_step(LammpcSession *session, const double *positions_ang,
                             const int *z, int n_atoms, const double *box9) {
  unsigned char step[16384];
  size_t step_size =
      make_step(positions_ang, z, n_atoms, box9, step, sizeof(step));
  size_t need = lammpc_potential_result_size_for_force_input(step, step_size);
  assert_true(need > 0);
  std::vector<unsigned char> out(need);
  size_t wrote = 0;
  LammpcResult r = lammpc_session_calculate_result(session, step, step_size,
                                                   out.data(), out.size(),
                                                   &wrote);
  return r;
}

static void assert_positions_bit_identical(LammpcSession *session,
                                           const double *expected,
                                           int n_atoms) {
  std::vector<double> got((size_t)n_atoms * 3u, 0.0);
  lammps_gather_atoms(session->lmp, "x", 1, 3, got.data());
  for (int i = 0; i < n_atoms * 3; ++i)
    assert_memory_equal(&got[i], &expected[i], sizeof(double));
}

static void test_upper_face_precision_orthogonal(void **state) {
  (void)state;
  LammpcSession *session = session_from_config("p p p");
  const double box9[9] = {10.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 10.0};
  const int z[2] = {14, 14};
  const double near = 10.0 - 1.0e-9;
  const double pos[6] = {0.125, 0.25, 0.5, near, near, near};

  LammpcResult r = run_step(session, pos, z, 2, box9);
  if (!r.ok)
    fail_msg("orthogonal near-face evaluation failed: %s", r.message);

  assert_int_equal((int)llround(lammps_get_natoms(session->lmp)), 2);
  assert_positions_bit_identical(session, pos, 2);
  lammpc_session_destroy(session);
}

static void test_upper_face_precision_triclinic(void **state) {
  (void)state;
  LammpcSession *session = session_from_config("p p p");
  const double box9[9] = {8.0, 0.0, 0.0, 1.0, 8.0, 0.0, 0.5, 0.25, 8.0};
  const int z[2] = {14, 14};
  const double eps = std::ldexp(1.0, -40);
  const double s = 1.0 - eps;
  const double pos[6] = {
      0.25, 0.5, 0.75,
      box9[0] * s + box9[3] * s + box9[6] * s,
      box9[4] * s + box9[7] * s,
      box9[8] * s,
  };

  LammpcResult r = run_step(session, pos, z, 2, box9);
  if (!r.ok)
    fail_msg("triclinic near-face evaluation failed: %s", r.message);

  assert_int_equal((int)llround(lammps_get_natoms(session->lmp)), 2);
  assert_positions_bit_identical(session, pos, 2);
  lammpc_session_destroy(session);
}

static void test_atom_count_mismatch_fails_loudly(void **state) {
  (void)state;
  LammpcSession *session = session_from_config("f f f");
  const double box9[9] = {10.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 10.0};
  const int z[2] = {14, 14};
  const double pos[6] = {0.1, 0.1, 0.1, 10.0, 5.0, 5.0};

  LammpcResult r = run_step(session, pos, z, 2, box9);
  assert_int_equal(r.ok, 0);
  assert_non_null(strstr(r.message, "atom count mismatch"));
  lammpc_session_destroy(session);
}

static int count_boundary_calls(void *user, int n_atoms, double *positions,
                                double *velocities,
                                const LammpcBoundaryBox *box) {
  (void)n_atoms;
  (void)positions;
  (void)velocities;
  (void)box;
  int *calls = static_cast<int *>(user);
  *calls += 1;
  return LAMMPC_BOUNDARY_ACCEPT;
}

static int shift_lj_boundary(void *user, int n_atoms, double *positions,
                             double *velocities,
                             const LammpcBoundaryBox *box) {
  (void)velocities;
  (void)box;
  int *calls = static_cast<int *>(user);
  *calls += 1;
  assert_int_equal(n_atoms, 2);
  positions[0] -= 0.05;
  positions[3] += 0.05;
  return LAMMPC_BOUNDARY_REFLECT;
}

static int reject_lj_boundary(void *user, int n_atoms, double *positions,
                              double *velocities,
                              const LammpcBoundaryBox *box) {
  (void)positions;
  (void)velocities;
  (void)box;
  int *calls = static_cast<int *>(user);
  *calls += 1;
  assert_int_equal(n_atoms, 2);
  return LAMMPC_BOUNDARY_REJECT;
}

static double lj_force(double separation) {
  const double sr6 = pow(2.5 / separation, 6.0);
  return 24.0 * 0.2 * (2.0 * sr6 * sr6 - sr6) / separation;
}

static void test_reflect_refreshes_lammps_forces(void **state) {
  (void)state;
  unsigned char config[16384];
  const size_t config_size = make_lj_config(config, sizeof(config));
  LammpcSession *session =
      lammpc_session_create_from_config(config, config_size);
  assert_non_null(session);

  const double box9[9] = {20.0, 0.0, 0.0, 0.0, 20.0,
                          0.0,  0.0, 0.0, 20.0};
  const int z[2] = {29, 29};
  const double initial_positions[6] = {5.0, 5.0, 5.0, 8.0, 5.0, 5.0};
  LammpcResult result =
      run_step(session, initial_positions, z, 2, box9);
  if (!result.ok)
    fail_msg("LJ reflection setup failed: %s", result.message);

  const double velocities[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  assert_int_equal(lammpc_session_nve_set_velocities(session, velocities, 6),
                   0);
  if (run_cmd(session->lmp, "timestep 0.000001") != 0 ||
      run_cmd(session->lmp, "fix 1 all nve") != 0)
    fail_msg("LJ reflection integrator setup failed: %s", lammpc_last_error());

  int calls = 0;
  assert_int_equal(
      lammpc_session_set_hard_boundary(session, shift_lj_boundary, &calls), 0);
  assert_int_equal(lammpc_session_integrate(session, 1), 0);
  assert_int_equal(calls, 1);

  double positions[6] = {};
  double forces[6] = {};
  lammps_gather_atoms(session->lmp, "x", 1, 3, positions);
  lammps_gather_atoms(session->lmp, "f", 1, 3, forces);
  const double separation = positions[3] - positions[0];
  const double expected_force = lj_force(separation);
  assert_true(fabs(forces[3] - expected_force) <
              fabs(expected_force) * 1.0e-10);
  assert_true(fabs(forces[0] + expected_force) <
              fabs(expected_force) * 1.0e-10);
  lammpc_session_destroy(session);
}

static void test_reject_refreshes_lammps_forces(void **state) {
  (void)state;
  unsigned char config[16384];
  const size_t config_size = make_lj_config(config, sizeof(config));
  LammpcSession *session =
      lammpc_session_create_from_config(config, config_size);
  assert_non_null(session);

  const double box9[9] = {20.0, 0.0, 0.0, 0.0, 20.0,
                          0.0,  0.0, 0.0, 20.0};
  const int z[2] = {29, 29};
  const double initial_positions[6] = {5.0, 5.0, 5.0, 8.0, 5.0, 5.0};
  LammpcResult result =
      run_step(session, initial_positions, z, 2, box9);
  if (!result.ok)
    fail_msg("LJ rejection setup failed: %s", result.message);

  const double velocities[6] = {-10.0, 0.0, 0.0, 10.0, 0.0, 0.0};
  assert_int_equal(lammpc_session_nve_set_velocities(session, velocities, 6),
                   0);
  if (run_cmd(session->lmp, "timestep 0.01") != 0 ||
      run_cmd(session->lmp, "fix 1 all nve") != 0 ||
      run_cmd(session->lmp, "compute test_pe_atom all pe/atom") != 0 ||
      run_cmd(session->lmp,
              "compute test_pe_sum all reduce sum c_test_pe_atom") != 0 ||
      run_cmd(session->lmp,
              "compute test_stress_atom all stress/atom NULL") != 0 ||
      run_cmd(session->lmp,
              "compute test_stress_sum all reduce sum "
              "c_test_stress_atom[1]") != 0 ||
      run_cmd(session->lmp,
              "thermo_style custom step pe pxx c_test_pe_sum "
              "c_test_stress_sum") != 0 ||
      run_cmd(session->lmp, "thermo 1") != 0)
    fail_msg("LJ rejection integrator setup failed: %s", lammpc_last_error());

  int calls = 0;
  assert_int_equal(
      lammpc_session_set_hard_boundary(session, reject_lj_boundary, &calls), 0);
  assert_int_equal(lammpc_session_integrate(session, 1), 0);
  assert_int_equal(calls, 1);

  double restored_positions[6] = {};
  double restored_forces[6] = {};
  lammps_gather_atoms(session->lmp, "x", 1, 3, restored_positions);
  lammps_gather_atoms(session->lmp, "f", 1, 3, restored_forces);
  assert_memory_equal(restored_positions, initial_positions,
                      sizeof(restored_positions));
  const double restored_pe = lammps_get_thermo(session->lmp, "pe");
  const double restored_pxx = lammps_get_thermo(session->lmp, "pxx");
  const double restored_pe_atom_sum =
      lammps_get_thermo(session->lmp, "c_test_pe_sum");
  const double restored_stress_atom_sum =
      lammps_get_thermo(session->lmp, "c_test_stress_sum");

  if (run_cmd(session->lmp, "run 0 pre yes post no") != 0)
    fail_msg("fresh force evaluation failed: %s", lammpc_last_error());
  double fresh_forces[6] = {};
  lammps_gather_atoms(session->lmp, "f", 1, 3, fresh_forces);
  constexpr double relative_tolerance = 1.0e-12;
  for (size_t i = 0; i < 6; ++i) {
    const double scale = std::max(1.0, fabs(fresh_forces[i]));
    assert_true(fabs(restored_forces[i] - fresh_forces[i]) <=
                relative_tolerance * scale);
  }
  assert_true(fabs(restored_pe - lammps_get_thermo(session->lmp, "pe")) <=
              relative_tolerance * std::max(1.0, fabs(restored_pe)));
  assert_true(fabs(restored_pxx - lammps_get_thermo(session->lmp, "pxx")) <=
              relative_tolerance * std::max(1.0, fabs(restored_pxx)));
  assert_true(fabs(restored_pe_atom_sum -
                   lammps_get_thermo(session->lmp, "c_test_pe_sum")) <=
              relative_tolerance * std::max(1.0, fabs(restored_pe_atom_sum)));
  assert_true(fabs(restored_stress_atom_sum -
                   lammps_get_thermo(session->lmp, "c_test_stress_sum")) <=
              relative_tolerance *
                  std::max(1.0, fabs(restored_stress_atom_sum)));
  lammpc_session_destroy(session);
}

/*
 * Destroy drops the host callback before it closes LAMMPS. Attach the
 * fix, prove one step calls the host, disarm, and prove the next step
 * does not. The fix stays installed across that second step, so a disarm
 * that does not clear the binding still invokes the host.
 */
static void test_destroy_disarms_hard_boundary(void **state) {
  (void)state;
  lammpc_embed_boundary_disarm(nullptr);

  LammpcSession *session = session_from_config("f f f");
  const double box9[9] = {10.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 10.0};
  const int z[2] = {14, 14};
  const double pos[6] = {1.0, 1.0, 1.0, 2.0, 2.0, 2.0};
  LammpcResult r = run_step(session, pos, z, 2, box9);
  if (!r.ok)
    fail_msg("disarm setup evaluation failed: %s", r.message);

  int calls = 0;
  assert_int_equal(
      lammpc_session_set_hard_boundary(session, count_boundary_calls, &calls),
      0);
  assert_int_equal(
      lammpc_embed_boundary_attach(session->lmp, count_boundary_calls, &calls),
      0);
  assert_int_equal(calls, 0);
  if (run_cmd(session->lmp, "timestep 0.001") != 0 ||
      run_cmd(session->lmp, "fix 1 all nve") != 0 ||
      run_cmd(session->lmp, "run 1") != 0)
    fail_msg("armed step failed: %s", lammpc_last_error());
  assert_int_equal(calls, 1);

  calls = 0;
  forget_hard_boundary(session);
  if (run_cmd(session->lmp, "run 1") != 0)
    fail_msg("disarmed step failed: %s", lammpc_last_error());
  assert_int_equal(calls, 0);
  lammpc_session_destroy(session);
}

int main(void) {
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_upper_face_precision_orthogonal),
      cmocka_unit_test(test_upper_face_precision_triclinic),
      cmocka_unit_test(test_atom_count_mismatch_fails_loudly),
      cmocka_unit_test(test_reflect_refreshes_lammps_forces),
      cmocka_unit_test(test_reject_refreshes_lammps_forces),
      cmocka_unit_test(test_destroy_disarms_hard_boundary),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
