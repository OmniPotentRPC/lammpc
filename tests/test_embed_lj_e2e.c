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

static size_t make_lj_config(const char *boundary, unsigned char *dst,
                             size_t dst_size) {
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

static size_t make_config(unsigned char *dst, size_t dst_size) {
  return make_lj_config("f f f", dst, dst_size);
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

static size_t make_boxed_step(const double coords[6], const double box9[9],
                              unsigned char *dst, size_t dst_size) {
  struct capn arena;
  capn_init_malloc(&arena);
  struct capn_segment *seg = capn_root(&arena).seg;

  ForceInput_ptr input = new_ForceInput(seg);
  capn_list64 pos = capn_new_list64(seg, 6);
  for (int i = 0; i < 6; ++i)
    capn_set64(pos, i, capn_from_f64(coords[i]));
  capn_list32 atmnrs = capn_new_list32(seg, 2);
  capn_set32(atmnrs, 0, 29);
  capn_set32(atmnrs, 1, 29);
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

static void evaluate_boxed_step(LammpcSession *session, const double coords[6],
                                const double box9[9], const char *context,
                                double *energy, double forces[6]) {
  unsigned char step[4096];
  size_t step_size = make_boxed_step(coords, box9, step, sizeof(step));
  size_t need = lammpc_potential_result_size_for_force_input(step, step_size);
  assert_true(need > 0);
  unsigned char *out = (unsigned char *)malloc(need);
  assert_non_null(out);
  size_t wrote = 0;
  LammpcResult r = lammpc_session_calculate_result(session, step, step_size,
                                                   out, need, &wrote);
  if (!r.ok)
    fail_msg("%s failed: %s", context, r.message);
  decode_result(out, wrote, energy, forces);
  free(out);
}

static void assert_energy_forces_match(const char *context, double e_ref,
                                       const double f_ref[6], double e_img,
                                       const double f_img[6]) {
  assert_true(fabs(e_ref - e_img) < 1e-12);
  for (int i = 0; i < 6; ++i) {
    if (fabs(f_ref[i] - f_img[i]) >= 1e-12)
      fail_msg("%s force mismatch at %d: ref=%.17g image=%.17g", context, i,
               f_ref[i], f_img[i]);
  }
}

static void fractional_to_cart(const double box9[9], const double s[3],
                               double r[3]) {
  r[0] = box9[0] * s[0] + box9[3] * s[1] + box9[6] * s[2];
  r[1] = box9[4] * s[1] + box9[7] * s[2];
  r[2] = box9[8] * s[2];
}

static void set_two_atom_positions(const double box9[9], const double s0[3],
                                   const double s1[3], double coords[6]) {
  double r0[3] = {0.0, 0.0, 0.0};
  double r1[3] = {0.0, 0.0, 0.0};
  fractional_to_cart(box9, s0, r0);
  fractional_to_cart(box9, s1, r1);
  coords[0] = r0[0];
  coords[1] = r0[1];
  coords[2] = r0[2];
  coords[3] = r1[0];
  coords[4] = r1[1];
  coords[5] = r1[2];
}

static void run_reused_session_image_case(const char *boundary,
                                          const double box9[9],
                                          const double reference[6],
                                          const double image[6],
                                          const char *name) {
  unsigned char config[8192];
  size_t config_size = make_lj_config(boundary, config, sizeof(config));
  LammpcSession *session =
      lammpc_session_create_from_config(config, config_size);
  assert_non_null(session);
  double e_ref = 0.0;
  double e_img = 0.0;
  double f_ref[6] = {0.0};
  double f_img[6] = {0.0};
  evaluate_boxed_step(session, reference, box9, "reference step", &e_ref, f_ref);
  evaluate_boxed_step(session, image, box9, name, &e_img, f_img);
  assert_energy_forces_match(name, e_ref, f_ref, e_img, f_img);
  lammpc_session_destroy(session);
}

static void run_creation_image_case(const char *boundary, const double box9[9],
                                    const double reference[6],
                                    const double image[6], const char *name) {
  unsigned char config[8192];
  size_t config_size = make_lj_config(boundary, config, sizeof(config));

  LammpcSession *session_ref =
      lammpc_session_create_from_config(config, config_size);
  LammpcSession *session_img =
      lammpc_session_create_from_config(config, config_size);
  assert_non_null(session_ref);
  assert_non_null(session_img);

  double e_ref = 0.0;
  double e_img = 0.0;
  double f_ref[6] = {0.0};
  double f_img[6] = {0.0};
  evaluate_boxed_step(session_ref, reference, box9, "creation reference", &e_ref,
                      f_ref);
  evaluate_boxed_step(session_img, image, box9, name, &e_img, f_img);
  assert_energy_forces_match(name, e_ref, f_ref, e_img, f_img);

  lammpc_session_destroy(session_ref);
  lammpc_session_destroy(session_img);
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

static void test_triclinic_plus2a_face_reproducer(void **state) {
  (void)state;
  const double box9[9] = {8.0, 0.0, 0.0, 1.0, 8.0, 0.0, 0.5, 0.25, 8.0};
  const double s1[3] = {0.45, 0.2, 0.35};
  const double s0_ref[3] = {0.0, 0.0, 0.3};
  const double s0_img[3] = {2.0, 0.0, 0.3};
  double ref[6] = {0.0};
  double img[6] = {0.0};
  set_two_atom_positions(box9, s0_ref, s1, ref);
  set_two_atom_positions(box9, s0_img, s1, img);
  run_reused_session_image_case("p p p", box9, ref, img,
                                "triclinic +2a face image");
}

static void test_triclinic_minus1a_face_reproducer(void **state) {
  (void)state;
  const double box9[9] = {10.86, 0.0, 0.0, 0.6, 10.86, 0.0, 0.3, 0.2, 10.86};
  const double s1[3] = {0.25, 0.35, 0.4};
  const double s0_ref[3] = {0.0, 1.0, 0.1};
  const double s0_img[3] = {-1.0, 1.0, 0.1};
  double ref[6] = {0.0};
  double img[6] = {0.0};
  set_two_atom_positions(box9, s0_ref, s1, ref);
  set_two_atom_positions(box9, s0_img, s1, img);
  run_reused_session_image_case("p p p", box9, ref, img,
                                "triclinic -1a face image");
}

static void test_orthogonal_nextafter_2l_reproducer(void **state) {
  (void)state;
  const double L = 7.992;
  const double box9[9] = {L, 0.0, 0.0, 0.0, L, 0.0, 0.0, 0.0, L};
  const double x_image = nextafter(2.0 * L, 0.0);
  const double x_reference = x_image - L;
  const double ref[6] = {x_reference, 1.5, 0.75, 2.3, 2.5, 1.0};
  const double img[6] = {x_image, 1.5, 0.75, 2.3, 2.5, 1.0};
  run_reused_session_image_case("p p p", box9, ref, img,
                                "orthogonal nextafter(2L,0) image");
}

static void test_partial_periodic_creation_ppf(void **state) {
  (void)state;
  const double box9[9] = {12.0, 0.0, 0.0, 1.5, 11.5, 0.0, 0.75, 0.5, 11.0};
  const double s0_ref[3] = {0.15, 0.2, 0.35};
  const double s1_ref[3] = {0.7, 0.6, 0.55};
  const double s0_img[3] = {3.15, -1.8, 0.35};
  const double s1_img[3] = {-3.3, 5.6, 0.55};
  double ref[6] = {0.0};
  double img[6] = {0.0};
  set_two_atom_positions(box9, s0_ref, s1_ref, ref);
  set_two_atom_positions(box9, s0_img, s1_img, img);
  run_creation_image_case("p p f", box9, ref, img, "creation p p f image");
}

static void test_partial_periodic_creation_pfp(void **state) {
  (void)state;
  const double box9[9] = {12.0, 0.0, 0.0, 1.5, 11.5, 0.0, 0.75, 0.5, 11.0};
  const double s0_ref[3] = {0.2, 0.35, 0.45};
  const double s1_ref[3] = {0.6, 0.55, 0.2};
  const double s0_img[3] = {4.2, 0.35, -2.55};
  const double s1_img[3] = {-3.4, 0.55, 5.2};
  double ref[6] = {0.0};
  double img[6] = {0.0};
  set_two_atom_positions(box9, s0_ref, s1_ref, ref);
  set_two_atom_positions(box9, s0_img, s1_img, img);
  run_creation_image_case("p f p", box9, ref, img, "creation p f p image");
}

static void test_partial_periodic_creation_fpp(void **state) {
  (void)state;
  const double box9[9] = {12.0, 0.0, 0.0, 1.5, 11.5, 0.0, 0.75, 0.5, 11.0};
  const double s0_ref[3] = {0.35, 0.25, 0.4};
  const double s1_ref[3] = {0.45, 0.7, 0.2};
  const double s0_img[3] = {0.35, -2.75, 4.4};
  const double s1_img[3] = {0.45, 5.7, -3.8};
  double ref[6] = {0.0};
  double img[6] = {0.0};
  set_two_atom_positions(box9, s0_ref, s1_ref, ref);
  set_two_atom_positions(box9, s0_img, s1_img, img);
  run_creation_image_case("f p p", box9, ref, img, "creation f p p image");
}

static void test_mixed_boundary_reuse_keeps_nonperiodic_axes(void **state) {
  (void)state;
  const double box9[9] = {12.0, 0.0, 0.0, 1.5, 11.5, 0.0, 0.75, 0.5, 11.0};
  struct boundary_case {
    const char *boundary;
    int periodic[3];
    int nonperiodic_axis;
  };
  const struct boundary_case cases[] = {
      {"p p s", {1, 1, 0}, 2},
      {"p s p", {1, 0, 1}, 1},
      {"s p p", {0, 1, 1}, 0},
  };

  for (size_t c = 0; c < sizeof(cases) / sizeof(cases[0]); ++c) {
    const struct boundary_case *item = &cases[c];
    double s0_ref[3] = {0.2, 0.3, 0.4};
    double s1_ref[3] = {0.7, 0.55, 0.25};
    s0_ref[item->nonperiodic_axis] = 1.25;
    s1_ref[item->nonperiodic_axis] = -0.25;
    double s0_img[3] = {s0_ref[0], s0_ref[1], s0_ref[2]};
    double s1_img[3] = {s1_ref[0], s1_ref[1], s1_ref[2]};
    const int shift0[3] = {4, -3, 5};
    const int shift1[3] = {-5, 6, -4};
    for (int d = 0; d < 3; ++d) {
      if (!item->periodic[d])
        continue;
      s0_img[d] += (double)shift0[d];
      s1_img[d] += (double)shift1[d];
    }
    double ref[6] = {0.0};
    double img[6] = {0.0};
    set_two_atom_positions(box9, s0_ref, s1_ref, ref);
    set_two_atom_positions(box9, s0_img, s1_img, img);
    run_reused_session_image_case(item->boundary, box9, ref, img,
                                  "mixed-boundary image");
  }
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
      cmocka_unit_test(test_triclinic_plus2a_face_reproducer),
      cmocka_unit_test(test_triclinic_minus1a_face_reproducer),
      cmocka_unit_test(test_orthogonal_nextafter_2l_reproducer),
      cmocka_unit_test(test_partial_periodic_creation_ppf),
      cmocka_unit_test(test_partial_periodic_creation_pfp),
      cmocka_unit_test(test_partial_periodic_creation_fpp),
      cmocka_unit_test(test_mixed_boundary_reuse_keeps_nonperiodic_axes),
      cmocka_unit_test(test_projected_md_stays_on_plane),
      cmocka_unit_test(test_overlay_rejected),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
