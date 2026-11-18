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

static size_t make_lj_config_with_units(const char *units, const char *boundary,
                                        unsigned char *dst, size_t dst_size) {
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
  pview.unitsStyle =
      (capn_text){(int)strlen(units), (char *)units, NULL};
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

static size_t make_lj_config(const char *boundary, unsigned char *dst,
                             size_t dst_size) {
  return make_lj_config_with_units("metal", boundary, dst, dst_size);
}

static size_t make_config(unsigned char *dst, size_t dst_size) {
  return make_lj_config("f f f", dst, dst_size);
}

static size_t make_step_positions_units(const double coordinates[6],
                                        const char *length_unit,
                                        const char *energy_unit,
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
  view.lengthUnit =
      (capn_text){(int)strlen(length_unit), (char *)length_unit, NULL};
  view.energyUnit =
      (capn_text){(int)strlen(energy_unit), (char *)energy_unit, NULL};
  write_ForceInput(&view, input);

  size_t written = write_flat(&arena, input.p, dst, dst_size);
  capn_free(&arena);
  return written;
}

static size_t make_step_positions(const double coordinates[6],
                                  unsigned char *dst, size_t dst_size) {
  return make_step_positions_units(coordinates, "angstrom", "eV", dst,
                                   dst_size);
}

static size_t make_step(double separation_ang, unsigned char *dst,
                        size_t dst_size) {
  const double coordinates[6] = {0.0, 0.0, 0.0,
                                 separation_ang, 0.0, 0.0};
  return make_step_positions(coordinates, dst, dst_size);
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

static void test_size_query_rejects_unsupported_energy_unit(void **state) {
  (void)state;
  const double coordinates[6] = {0.0, 0.0, 0.0, 2.5, 0.0, 0.0};
  unsigned char step[4096];
  size_t good = make_step_positions_units(coordinates, "angstrom", "eV", step,
                                sizeof(step));
  assert_true(lammpc_potential_result_size_for_force_input(step, good) > 0);
  size_t bad =
      make_step_positions_units(coordinates, "angstrom", "furlong", step,
                                sizeof(step));
  assert_true(bad > 0);
  assert_int_equal(lammpc_potential_result_size_for_force_input(step, bad), 0);
}

static const size_t kSizeSentinel = 99;

static void assert_failed_with(LammpcResult result, const char *needle) {
  assert_int_equal(result.ok, 0);
  assert_non_null(strstr(result.message, needle));
  assert_non_null(strstr(lammpc_last_error(), needle));
}

static void test_session_size_zero_for_unsupported_units(void **state) {
  (void)state;
  unsigned char config[8192];
  size_t config_size = make_config(config, sizeof(config));
  LammpcSession *session =
      lammpc_session_create_from_config(config, config_size);
  assert_non_null(session);

  const double coordinates[6] = {0.0, 0.0, 0.0, 2.5, 0.0, 0.0};
  const char *const cases[2][2] = {{"furlong", "eV"}, {"angstrom", "furlong"}};
  for (size_t i = 0; i < 2; ++i) {
    unsigned char step[4096];
    size_t step_size = make_step_positions_units(
        coordinates, cases[i][0], cases[i][1], step, sizeof(step));
    unsigned char out[4096];
    size_t size = kSizeSentinel;
    LammpcResult r = lammpc_session_calculate_result(
        session, step, step_size, out, sizeof(out), &size);
    assert_failed_with(r, "invalid ForceInput geometry or units");
    assert_int_equal(size, 0);

    size = kSizeSentinel;
    r = lammpc_session_calculate_result(session, step, step_size, NULL, 0,
                                        &size);
    assert_failed_with(r, "invalid ForceInput geometry or units");
    assert_int_equal(size, 0);
  }
  lammpc_session_destroy(session);
}

static void test_session_size_zero_for_unreadable_force_input(void **state) {
  (void)state;
  unsigned char config[8192];
  size_t config_size = make_config(config, sizeof(config));
  LammpcSession *session =
      lammpc_session_create_from_config(config, config_size);
  assert_non_null(session);

  const unsigned char junk[3] = {1, 2, 3};
  unsigned char out[4096];
  size_t size = kSizeSentinel;
  LammpcResult r = lammpc_session_calculate_result(
      session, junk, sizeof(junk), out, sizeof(out), &size);
  assert_failed_with(r, "invalid ForceInput message");
  assert_int_equal(size, 0);
  lammpc_session_destroy(session);
}

static void test_session_size_reports_required_on_small_buffer(void **state) {
  (void)state;
  unsigned char config[8192];
  size_t config_size = make_config(config, sizeof(config));
  unsigned char step[4096];
  size_t step_size = make_step(3.0, step, sizeof(step));
  LammpcSession *session =
      lammpc_session_create_from_config(config, config_size);
  assert_non_null(session);

  size_t need = lammpc_potential_result_size_for_force_input(step, step_size);
  assert_true(need > 0);
  unsigned char out[4096];
  assert_true(need > 1 && need <= sizeof(out));

  size_t size = kSizeSentinel;
  LammpcResult r = lammpc_session_calculate_result(session, step, step_size,
                                                   NULL, 0, &size);
  assert_failed_with(r, "PotentialResult buffer too small");
  assert_int_equal(size, need);

  size = kSizeSentinel;
  r = lammpc_session_calculate_result(session, step, step_size, out, need - 1,
                                      &size);
  assert_failed_with(r, "PotentialResult buffer too small");
  assert_int_equal(size, need);

  size = kSizeSentinel;
  r = lammpc_session_calculate_result(session, step, step_size, out, need,
                                      &size);
  assert_int_equal(r.ok, 1);
  assert_true(size > 0 && size <= need);
  lammpc_session_destroy(session);
}

static void test_session_null_size_pointer_is_invalid_arguments(void **state) {
  (void)state;
  unsigned char config[8192];
  size_t config_size = make_config(config, sizeof(config));
  unsigned char step[4096];
  size_t step_size = make_step(3.0, step, sizeof(step));
  LammpcSession *session =
      lammpc_session_create_from_config(config, config_size);
  assert_non_null(session);

  unsigned char out[4096];
  LammpcResult r = lammpc_session_calculate_result(
      session, step, step_size, out, sizeof(out), NULL);
  assert_failed_with(r, "invalid arguments");
  lammpc_session_destroy(session);
}

static void test_one_shot_size_zero_for_unsupported_units(void **state) {
  (void)state;
  unsigned char config[8192];
  size_t config_size = make_config(config, sizeof(config));
  const double coordinates[6] = {0.0, 0.0, 0.0, 2.5, 0.0, 0.0};
  unsigned char step[4096];
  size_t step_size = make_step_positions_units(coordinates, "angstrom",
                                               "furlong", step, sizeof(step));
  unsigned char out[4096];
  size_t size = kSizeSentinel;
  LammpcResult r = lammpc_calculate_result_from_config(
      config, config_size, step, step_size, out, sizeof(out), &size);
  assert_failed_with(r, "invalid ForceInput geometry or units");
  assert_int_equal(size, 0);
}

static void test_one_shot_config_creation_failure_zeroes_size(void **state) {
  (void)state;
  unsigned char bad_config[8192];
  size_t bad_size =
      make_lj_config_with_units("real", "f f f", bad_config, sizeof(bad_config));
  unsigned char step[4096];
  size_t step_size = make_step(3.0, step, sizeof(step));
  unsigned char out[4096];

  size_t size = kSizeSentinel;
  LammpcResult r = lammpc_calculate_result_from_config(
      bad_config, bad_size, step, step_size, out, sizeof(out), &size);
  assert_failed_with(r, "only unitsStyle=metal is supported");
  assert_int_equal(size, 0);

  r = lammpc_calculate_result_from_config(bad_config, bad_size, step,
                                          step_size, out, sizeof(out), NULL);
  assert_failed_with(r, "invalid arguments");
  assert_null(strstr(lammpc_last_error(), "unitsStyle"));
}

static void test_one_shot_params_creation_failure_zeroes_size(void **state) {
  (void)state;
  unsigned char step[4096];
  size_t step_size = make_step(3.0, step, sizeof(step));
  unsigned char out[4096];

  size_t size = kSizeSentinel;
  LammpcResult r = lammpc_calculate_result(NULL, 0, step, step_size, out,
                                           sizeof(out), &size);
  assert_failed_with(r, "invalid LammpsParams message");
  assert_int_equal(size, 0);

  r = lammpc_calculate_result(NULL, 0, step, step_size, out, sizeof(out),
                              NULL);
  assert_failed_with(r, "invalid arguments");
  assert_null(strstr(lammpc_last_error(), "LammpsParams"));
}

static void test_one_shot_params_size_zero_for_unsupported_units(void **state) {
  (void)state;
  struct capn arena;
  capn_init_malloc(&arena);
  struct capn_segment *seg = capn_root(&arena).seg;
  LammpsParams_ptr params = new_LammpsParams(seg);
  capn_ptr coeffs = capn_new_ptr_list(seg, 1);
  assert_int_equal(
      capn_set_text(coeffs, 0, (capn_text){11, "1 1 0.2 2.5", NULL}), 0);
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
  unsigned char params_bytes[8192];
  size_t params_size = write_flat(&arena, params.p, params_bytes,
                                  sizeof(params_bytes));
  capn_free(&arena);

  const double coordinates[6] = {0.0, 0.0, 0.0, 2.5, 0.0, 0.0};
  unsigned char step[4096];
  size_t step_size = make_step_positions_units(coordinates, "furlong", "eV",
                                               step, sizeof(step));
  unsigned char out[4096];
  size_t size = kSizeSentinel;
  LammpcResult r = lammpc_calculate_result(params_bytes, params_size, step,
                                           step_size, out, sizeof(out), &size);
  assert_failed_with(r, "invalid ForceInput geometry or units");
  assert_int_equal(size, 0);
}

static void test_real_units_are_rejected(void **state) {
  (void)state;
  unsigned char config[8192];
  const size_t config_size =
      make_lj_config_with_units("real", "f f f", config, sizeof(config));

  assert_null(lammpc_session_create_from_config(config, config_size));
  assert_string_equal(lammpc_last_error(),
                      "only unitsStyle=metal is supported");
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

static void test_md_two_session_replay_is_byte_identical(void **state) {
  (void)state;
  unsigned char config[8192];
  const size_t config_size = make_config(config, sizeof(config));
  unsigned char step[4096];
  const size_t step_size = make_step(3.0, step, sizeof(step));
  LammpcSession *first =
      lammpc_session_create_from_config(config, config_size);
  LammpcSession *second =
      lammpc_session_create_from_config(config, config_size);
  assert_non_null(first);
  assert_non_null(second);

  const double inverse_sqrt_two = 1.0 / sqrt(2.0);
  const double point[6] = {0.0, 0.0, 0.0, 3.0, 0.0, 0.0};
  const double normal[6] = {-inverse_sqrt_two, 0.0, 0.0,
                            inverse_sqrt_two,  0.0, 0.0};
  const LammpcMdConfig controls = {
      .time_step_ps = 0.0002,
      .thermal_energy_ev = 0.01,
      .friction_per_ps = 10.0,
      .seed = 41,
      .coordinates = LAMMPC_PLANE_MASS_WEIGHTED,
  };
  assert_int_equal(lammpc_session_md_initialize(
                       first, step, step_size, point, normal, 6, &controls),
                   0);
  assert_int_equal(lammpc_session_md_initialize(
                       second, step, step_size, point, normal, 6, &controls),
                   0);
  double masses[2] = {0.0, 0.0};
  assert_int_equal(lammpc_session_masses(first, masses, 2), 0);
  assert_true(masses[0] > 0.0);
  assert_true(fabs(masses[0] - masses[1]) < 1.0e-12);
  const double initial_velocities[6] = {0.0, 0.25, 0.0, 0.0, -0.25, 0.0};
  assert_int_equal(lammpc_session_md_set_velocities(
                       first, initial_velocities, 6, controls.seed),
                   0);
  assert_int_equal(lammpc_session_md_set_velocities(
                       second, initial_velocities, 6, controls.seed),
                   0);
  assert_int_equal(lammpc_session_md_step(first, 64), 0);
  assert_int_equal(lammpc_session_md_step(second, 64), 0);

  double first_positions[6], first_velocities[6], first_forces[6];
  double first_removed[6], second_positions[6], second_velocities[6];
  double second_forces[6], second_removed[6];
  LammpcMdObservation first_observation, second_observation;
  assert_int_equal(lammpc_session_md_state(
                       first, first_positions, first_velocities, first_forces,
                       first_removed, 6, &first_observation),
                   0);
  assert_int_equal(lammpc_session_md_state(
                       second, second_positions, second_velocities,
                       second_forces, second_removed, 6, &second_observation),
                   0);
  assert_memory_equal(first_positions, second_positions,
                      sizeof(first_positions));
  assert_memory_equal(first_velocities, second_velocities,
                      sizeof(first_velocities));
  assert_memory_equal(first_forces, second_forces, sizeof(first_forces));
  assert_memory_equal(first_removed, second_removed, sizeof(first_removed));
  assert_memory_equal(&first_observation.constraint_force,
                      &second_observation.constraint_force, sizeof(double));
  assert_int_equal(first_observation.force_evaluations,
                   second_observation.force_evaluations);
  assert_int_equal(first_observation.steps, second_observation.steps);
  lammpc_session_destroy(first);
  lammpc_session_destroy(second);
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

static void test_early_returns_refresh_last_error(void **state) {
  (void)state;

  assert_int_not_equal(lammpc_set_params(NULL, 0), 0);
  assert_non_null(strstr(lammpc_last_error(),
                         "global params are not supported"));

  assert_null(lammpc_session_create_from_config(NULL, 0));
  assert_non_null(strstr(lammpc_last_error(),
                         "PotentialConfig bytes are required"));
  assert_null(strstr(lammpc_last_error(),
                     "global params are not supported"));

  assert_int_not_equal(lammpc_set_params(NULL, 0), 0);
  assert_non_null(strstr(lammpc_last_error(),
                         "global params are not supported"));
  assert_int_not_equal(lammpc_session_set_params(NULL, NULL, 0), 0);
  assert_non_null(strstr(lammpc_last_error(),
                         "session_set_params requires a session"));
  assert_null(strstr(lammpc_last_error(),
                     "global params are not supported"));

  assert_int_not_equal(lammpc_set_params(NULL, 0), 0);
  assert_non_null(strstr(lammpc_last_error(),
                         "global params are not supported"));
  assert_int_not_equal(lammpc_session_configure(NULL, NULL, 0), 0);
  assert_non_null(strstr(lammpc_last_error(),
                         "session_configure requires a session"));
  assert_null(strstr(lammpc_last_error(),
                     "global params are not supported"));
}

static void test_calculate_result_null_refreshes_last_error(void **state) {
  (void)state;

  assert_int_not_equal(lammpc_set_params(NULL, 0), 0);
  assert_non_null(strstr(lammpc_last_error(),
                         "global params are not supported"));

  size_t out_size = 99;
  LammpcResult result = lammpc_session_calculate_result(NULL, NULL, 0, NULL, 0,
                                                        &out_size);
  assert_int_equal(result.ok, 0);
  assert_non_null(strstr(result.message, "invalid arguments"));
  assert_non_null(strstr(lammpc_last_error(), "invalid arguments"));
  assert_null(strstr(lammpc_last_error(), "global params are not supported"));
  assert_int_equal(out_size, 99);

  assert_int_not_equal(lammpc_set_params(NULL, 0), 0);
  LammpcResult forces = lammpc_energy_forces(0, NULL, NULL, NULL, 0, NULL);
  assert_int_equal(forces.ok, 0);
  assert_non_null(strstr(lammpc_last_error(),
                         "raw-array calls need a session"));
  assert_null(strstr(lammpc_last_error(), "global params are not supported"));

  assert_int_not_equal(lammpc_set_params(NULL, 0), 0);
  LammpcResult energy = lammpc_energy(0, NULL, NULL, NULL, 0);
  assert_int_equal(energy.ok, 0);
  assert_non_null(strstr(lammpc_last_error(),
                         "raw-array calls need a session"));

  assert_int_not_equal(lammpc_set_params(NULL, 0), 0);
  LammpcResult gradient = lammpc_energy_gradient(0, NULL, NULL, NULL, 0, NULL);
  assert_int_equal(gradient.ok, 0);
  assert_non_null(strstr(lammpc_last_error(),
                         "raw-array calls need a session"));

  unsigned char config[8192];
  size_t config_size = make_config(config, sizeof(config));
  assert_int_not_equal(lammpc_set_params(NULL, 0), 0);
  LammpcResult one_shot = lammpc_calculate_result_from_config(
      config, config_size, NULL, 0, NULL, 0, NULL);
  assert_int_equal(one_shot.ok, 0);
  assert_non_null(strstr(one_shot.message, "invalid arguments"));
  assert_non_null(strstr(lammpc_last_error(), "invalid arguments"));
  assert_null(strstr(lammpc_last_error(), "global params are not supported"));

  unsigned char step[4096];
  size_t step_size = make_step(3.0, step, sizeof(step));
  LammpcSession *session =
      lammpc_session_create_from_config(config, config_size);
  assert_non_null(session);
  assert_int_not_equal(lammpc_set_params(NULL, 0), 0);
  size_t need = lammpc_potential_result_size_for_force_input(step, step_size);
  assert_true(need > 0);
  unsigned char *out = (unsigned char *)malloc(need);
  assert_non_null(out);
  size_t wrote = 0;
  LammpcResult ok = lammpc_session_calculate_result(session, step, step_size,
                                                    out, need, &wrote);
  free(out);
  lammpc_session_destroy(session);
  assert_int_equal(ok.ok, 1);
  assert_string_equal(ok.message, "ok");
  assert_string_equal(lammpc_last_error(), "");
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
      cmocka_unit_test(test_real_units_are_rejected),
      cmocka_unit_test(test_session_size_zero_for_unsupported_units),
      cmocka_unit_test(test_session_size_zero_for_unreadable_force_input),
      cmocka_unit_test(test_session_size_reports_required_on_small_buffer),
      cmocka_unit_test(test_session_null_size_pointer_is_invalid_arguments),
      cmocka_unit_test(test_one_shot_size_zero_for_unsupported_units),
      cmocka_unit_test(test_one_shot_config_creation_failure_zeroes_size),
      cmocka_unit_test(test_one_shot_params_creation_failure_zeroes_size),
      cmocka_unit_test(test_one_shot_params_size_zero_for_unsupported_units),
      cmocka_unit_test(test_size_query_rejects_unsupported_energy_unit),
      cmocka_unit_test(test_triclinic_plus2a_face_reproducer),
      cmocka_unit_test(test_triclinic_minus1a_face_reproducer),
      cmocka_unit_test(test_orthogonal_nextafter_2l_reproducer),
      cmocka_unit_test(test_partial_periodic_creation_ppf),
      cmocka_unit_test(test_partial_periodic_creation_pfp),
      cmocka_unit_test(test_partial_periodic_creation_fpp),
      cmocka_unit_test(test_mixed_boundary_reuse_keeps_nonperiodic_axes),
      cmocka_unit_test(test_projected_md_stays_on_plane),
      cmocka_unit_test(test_md_two_session_replay_is_byte_identical),
      cmocka_unit_test(test_md_step_applies_hard_boundary_each_step),
      cmocka_unit_test(test_md_reject_restores_pre_step_force_state),
      cmocka_unit_test(test_md_invalid_decision_restores_pre_step_state),
      cmocka_unit_test(test_md_reflect_reprojects_and_reevaluates),
      cmocka_unit_test(test_early_returns_refresh_last_error),
      cmocka_unit_test(test_calculate_result_null_refreshes_last_error),
      cmocka_unit_test(test_overlay_rejected),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
