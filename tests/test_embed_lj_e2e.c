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

static size_t make_step(double separation_ang, unsigned char *dst,
                        size_t dst_size) {
  struct capn arena;
  capn_init_malloc(&arena);
  struct capn_segment *seg = capn_root(&arena).seg;

  ForceInput_ptr input = new_ForceInput(seg);
  capn_list64 pos = capn_new_list64(seg, 6);
  const double coords[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  for (int i = 0; i < 6; ++i)
    capn_set64(pos, i, capn_from_f64(coords[i]));
  capn_set64(pos, 3, capn_from_f64(separation_ang));
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
      cmocka_unit_test(test_overlay_rejected),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
