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
#include <unistd.h>

#include <cmocka.h>

static size_t write_flat(struct capn *arena, capn_ptr root_struct,
                         unsigned char *dst, size_t dst_size) {
  capn_ptr root = capn_root(arena);
  assert_int_equal(capn_setp(root, 0, root_struct), 0);
  int written = capn_write_mem(arena, dst, dst_size, 0);
  assert_true(written > 0);
  return (size_t)written;
}

static size_t make_cu_eam_config(const char *potential_path, unsigned char *dst,
                                 size_t dst_size) {
  struct capn arena;
  capn_init_malloc(&arena);
  struct capn_segment *seg = capn_root(&arena).seg;

  char coeff_line[2048];
  snprintf(coeff_line, sizeof(coeff_line), "* * %s Cu", potential_path);

  LammpsParams_ptr params = new_LammpsParams(seg);
  capn_ptr coeffs = capn_new_ptr_list(seg, 1);
  capn_text coeff = {(int)strlen(coeff_line), coeff_line, 0};
  assert_int_equal(capn_set_text(coeffs, 0, coeff), 0);
  capn_list32 type_map = capn_new_list32(seg, 1);
  capn_set32(type_map, 0, 29);
  capn_list64 masses = capn_new_list64(seg, 1);
  capn_set64(masses, 0, capn_from_f64(63.546));

  struct LammpsParams pview;
  memset(&pview, 0, sizeof(pview));
  pview.unitsStyle = (capn_text){5, "metal", NULL};
  pview.atomStyle = (capn_text){6, "atomic", NULL};
  pview.pairStyle = (capn_text){9, "eam/alloy", NULL};
  pview.pairCoeffs = coeffs;
  pview.typeToAtomicNumber = type_map;
  pview.masses = masses;
  pview.newtonPair = 1;
  pview.boundary = (capn_text){5, "p p p", NULL};
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

static size_t make_force_input(const double *positions_ang, const int *z,
                               int n_atoms, const double *box9,
                               unsigned char *dst, size_t dst_size) {
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
  capn_list64 box = capn_new_list64(seg, 9);
  for (int i = 0; i < 9; ++i)
    capn_set64(box, i, capn_from_f64(box9[i]));

  struct ForceInput view;
  memset(&view, 0, sizeof(view));
  view.pos = pos;
  view.atmnrs = atmnrs;
  view.box = box;
  view.lengthUnit = (capn_text){8, "angstrom", NULL};
  view.energyUnit = (capn_text){2, "eV", NULL};
  write_ForceInput(&view, input);

  size_t written = write_flat(&arena, input.p, dst, dst_size);
  capn_free(&arena);
  return written;
}

static void build_fcc_2x2x2_cu(double *positions_ang, int *z, int n_atoms,
                               double *box9) {
  const double a = 3.615;
  const double basis[4][3] = {
      {0.0, 0.0, 0.0},
      {0.0, 0.5, 0.5},
      {0.5, 0.0, 0.5},
      {0.5, 0.5, 0.0},
  };
  int idx = 0;
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      for (int k = 0; k < 2; ++k) {
        for (int b = 0; b < 4; ++b) {
          positions_ang[3 * idx + 0] = (i + basis[b][0]) * a;
          positions_ang[3 * idx + 1] = (j + basis[b][1]) * a;
          positions_ang[3 * idx + 2] = (k + basis[b][2]) * a;
          z[idx] = 29;
          ++idx;
        }
      }
    }
  }
  assert_int_equal(idx, n_atoms);
  memset(box9, 0, sizeof(double) * 9u);
  box9[0] = 2.0 * a;
  box9[4] = 2.0 * a;
  box9[8] = 2.0 * a;
}

static void make_initial_velocities(double *velocities, int n_atoms) {
  double sum[3] = {0.0, 0.0, 0.0};
  for (int i = 0; i < n_atoms; ++i) {
    const double scale = 0.002 + 0.0002 * (double)(i % 5);
    velocities[3 * i + 0] = scale * (double)(i % 3 - 1);
    velocities[3 * i + 1] = scale * (double)((2 * i) % 3 - 1);
    velocities[3 * i + 2] = scale * (double)((3 * i + 1) % 3 - 1);
    sum[0] += velocities[3 * i + 0];
    sum[1] += velocities[3 * i + 1];
    sum[2] += velocities[3 * i + 2];
  }
  for (int i = 0; i < n_atoms; ++i) {
    velocities[3 * i + 0] -= sum[0] / (double)n_atoms;
    velocities[3 * i + 1] -= sum[1] / (double)n_atoms;
    velocities[3 * i + 2] -= sum[2] / (double)n_atoms;
  }
}

static void initialize_geometry(LammpcSession *session, const double *positions,
                                const int *z, int n_atoms,
                                const double *box9) {
  unsigned char step[32768];
  size_t step_size =
      make_force_input(positions, z, n_atoms, box9, step, sizeof(step));
  const size_t need =
      lammpc_potential_result_size_for_force_input(step, step_size);
  assert_true(need > 0);
  unsigned char *out = malloc(need);
  assert_non_null(out);
  size_t wrote = 0;
  LammpcResult result = lammpc_session_calculate_result(
      session, step, step_size, out, need, &wrote);
  if (!result.ok)
    fail_msg("session initialization failed: %s", result.message);
  free(out);
}

static double max_abs_delta(const double *a, const double *b, size_t n) {
  double max_delta = 0.0;
  for (size_t i = 0; i < n; ++i) {
    const double delta = fabs(a[i] - b[i]);
    if (delta > max_delta)
      max_delta = delta;
  }
  return max_delta;
}

static double max_abs_periodic_delta(const double *a, const double *b, int n_atoms,
                                     const double *box9) {
  const double box[3] = {box9[0], box9[4], box9[8]};
  double max_delta = 0.0;
  for (int i = 0; i < n_atoms; ++i) {
    for (int d = 0; d < 3; ++d) {
      const double raw = a[3 * i + d] - b[3 * i + d];
      const double wrapped = raw - nearbyint(raw / box[d]) * box[d];
      const double delta = fabs(wrapped);
      if (delta > max_delta)
        max_delta = delta;
    }
  }
  return max_delta;
}

static void assert_ok_with_error(int rc, const char *context) {
  if (rc != 0)
    fail_msg("%s: %s", context, lammpc_last_error());
}

static void test_nve_energy_and_time_reversibility(void **state) {
  (void)state;
  assert_int_equal(lammpc_available(), 1);
  const char *lammps_root = getenv("LAMMPS_ROOT");
  assert_non_null(lammps_root);
  char potential[1024];
  snprintf(potential, sizeof(potential),
           "%s/share/lammps/potentials/Cu_mishin1.eam.alloy", lammps_root);
  assert_int_equal(access(potential, R_OK), 0);

  enum { n_atoms = 32 };
  double positions0[n_atoms * 3];
  double box9[9];
  int z[n_atoms];
  build_fcc_2x2x2_cu(positions0, z, n_atoms, box9);
  double velocities0[n_atoms * 3];
  make_initial_velocities(velocities0, n_atoms);

  unsigned char config[32768];
  size_t config_size =
      make_cu_eam_config(potential, config, sizeof(config));
  LammpcSession *session =
      lammpc_session_create_from_config(config, config_size);
  assert_non_null(session);

  initialize_geometry(session, positions0, z, n_atoms, box9);
  assert_ok_with_error(lammpc_session_nve_set_positions(
                           session, positions0, n_atoms * 3u, box9, 9),
                       "nve_set_positions");
  assert_ok_with_error(
      lammpc_session_nve_set_velocities(session, velocities0, n_atoms * 3u),
      "nve_set_velocities");

  LammpcNveObservation obs0 = {0};
  assert_ok_with_error(
      lammpc_session_nve_state(session, NULL, NULL, n_atoms * 3u, &obs0),
      "nve_state(obs0)");
  const double total0 = obs0.potential_energy_ev + obs0.kinetic_energy_ev;

  const size_t n_steps = 300;
  const double dt_ps = 0.0005;
  assert_ok_with_error(lammpc_session_nve_step(session, n_steps, dt_ps),
                       "nve_step(forward)");

  double positions1[n_atoms * 3];
  double velocities1[n_atoms * 3];
  LammpcNveObservation obs1 = {0};
  assert_ok_with_error(
      lammpc_session_nve_state(session, positions1, velocities1, n_atoms * 3u,
                               &obs1),
      "nve_state(obs1)");
  const double total1 = obs1.potential_energy_ev + obs1.kinetic_energy_ev;
  assert_true(fabs(total1 - total0) < 2e-3);

  double reversed_velocities[n_atoms * 3];
  for (int i = 0; i < n_atoms * 3; ++i)
    reversed_velocities[i] = -velocities1[i];
  assert_ok_with_error(lammpc_session_nve_set_velocities(
                           session, reversed_velocities, n_atoms * 3u),
                       "nve_set_velocities(reverse)");
  assert_ok_with_error(lammpc_session_nve_step(session, n_steps, dt_ps),
                       "nve_step(reverse)");

  double positions2[n_atoms * 3];
  double velocities2[n_atoms * 3];
  LammpcNveObservation obs2 = {0};
  assert_ok_with_error(
      lammpc_session_nve_state(session, positions2, velocities2, n_atoms * 3u,
                               &obs2),
      "nve_state(obs2)");

  assert_true(max_abs_periodic_delta(positions2, positions0, n_atoms, box9) <
              3e-4);
  for (int i = 0; i < n_atoms * 3; ++i)
    reversed_velocities[i] = -velocities0[i];
  assert_true(max_abs_delta(velocities2, reversed_velocities, n_atoms * 3u) <
              3e-4);

  lammpc_session_destroy(session);
}

int main(void) {
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_nve_energy_and_time_reversibility),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
