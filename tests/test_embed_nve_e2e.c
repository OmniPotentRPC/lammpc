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

static size_t make_potential_config(const char *pair_style,
                                    const char *coeff_line,
                                    const char *boundary, unsigned char *dst,
                                    size_t dst_size) {
  struct capn arena;
  capn_init_malloc(&arena);
  struct capn_segment *seg = capn_root(&arena).seg;

  LammpsParams_ptr params = new_LammpsParams(seg);
  capn_ptr coeffs = capn_new_ptr_list(seg, 1);
  capn_text coeff = {(int)strlen(coeff_line), (char *)coeff_line, 0};
  assert_int_equal(capn_set_text(coeffs, 0, coeff), 0);
  capn_list32 type_map = capn_new_list32(seg, 1);
  capn_set32(type_map, 0, 29);
  capn_list64 masses = capn_new_list64(seg, 1);
  capn_set64(masses, 0, capn_from_f64(63.546));

  struct LammpsParams pview;
  memset(&pview, 0, sizeof(pview));
  pview.unitsStyle = (capn_text){5, "metal", NULL};
  pview.atomStyle = (capn_text){6, "atomic", NULL};
  pview.pairStyle =
      (capn_text){(int)strlen(pair_style), (char *)pair_style, NULL};
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

static size_t make_cu_eam_config(const char *potential_path, unsigned char *dst,
                                 size_t dst_size) {
  char coeff_line[2048];
  snprintf(coeff_line, sizeof(coeff_line), "* * %s Cu", potential_path);
  return make_potential_config("eam/alloy", coeff_line, "p p p", dst, dst_size);
}

static size_t make_lj_config(const char *boundary, unsigned char *dst,
                             size_t dst_size) {
  return make_potential_config("lj/cut 10.0", "1 1 0.2 2.5", boundary, dst,
                               dst_size);
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

static const double kMassVelocityToEv = 1.0364269656262175e-4;
static const double kCuMassAmu = 63.546;

static void fractional_to_cart(const double box9[9], const double s[3],
                               double r[3]) {
  r[0] = box9[0] * s[0] + box9[3] * s[1] + box9[6] * s[2];
  r[1] = box9[4] * s[1] + box9[7] * s[2];
  r[2] = box9[8] * s[2];
}

static void set_two_atom_positions(const double box9[9], const double s0[3],
                                   const double s1[3], double coords[6]) {
  fractional_to_cart(box9, s0, &coords[0]);
  fractional_to_cart(box9, s1, &coords[3]);
}

static void cart_to_fractional(const double box9[9], const double r[3],
                               double s[3]) {
  const double xprd = box9[0];
  const double yprd = box9[4];
  const double zprd = box9[8];
  const double xy = box9[3];
  const double xz = box9[6];
  const double yz = box9[7];
  s[0] = (1.0 / xprd) * r[0] + (-xy / (xprd * yprd)) * r[1] +
         ((yz * xy - yprd * xz) / (xprd * yprd * zprd)) * r[2];
  s[1] = (1.0 / yprd) * r[1] + (-yz / (yprd * zprd)) * r[2];
  s[2] = (1.0 / zprd) * r[2];
}

static void assert_positions_folded_and_matched(
    const char *boundary, const double box9[9], const double pos_ref[6],
    const double pos_img[6], const char *name) {
  char f0[8] = {0}, f1[8] = {0}, f2[8] = {0};
  assert_int_equal(sscanf(boundary, "%7s %7s %7s", f0, f1, f2), 3);
  const int periodic[3] = {f0[0] == 'p', f1[0] == 'p', f2[0] == 'p'};
  for (int a = 0; a < 2; ++a) {
    double s_ref[3] = {0.0}, s_img[3] = {0.0};
    cart_to_fractional(box9, &pos_ref[3 * a], s_ref);
    cart_to_fractional(box9, &pos_img[3 * a], s_img);
    for (int d = 0; d < 3; ++d) {
      if (periodic[d]) {
        assert_true(s_img[d] >= -1e-12 && s_img[d] <= 1.0 + 1e-12);
        const double ds = s_img[d] - s_ref[d];
        assert_true(fabs(ds - nearbyint(ds)) < 1e-12);
      } else {
        assert_true(fabs(s_img[d] - s_ref[d]) < 1e-12);
      }
    }
  }
  (void)name;
}

static void run_lj_nve_set_positions_case(const char *boundary,
                                          const double box9[9],
                                          const double reference[6],
                                          const double image[6],
                                          const char *name) {
  unsigned char config[8192];
  const size_t config_size = make_lj_config(boundary, config, sizeof(config));
  const int z[2] = {29, 29};
  const double s0_seed[3] = {0.25, 0.25, 0.25};
  const double s1_seed[3] = {0.55, 0.55, 0.55};
  double seed[6] = {0.0};
  set_two_atom_positions(box9, s0_seed, s1_seed, seed);

  for (int pass_box = 0; pass_box <= 1; ++pass_box) {
    LammpcSession *session_ref =
        lammpc_session_create_from_config(config, config_size);
    LammpcSession *session_img =
        lammpc_session_create_from_config(config, config_size);
    assert_non_null(session_ref);
    assert_non_null(session_img);

    initialize_geometry(session_ref, seed, z, 2, box9);
    initialize_geometry(session_img, seed, z, 2, box9);

    const double *setter_box = pass_box ? box9 : NULL;
    const size_t setter_box_len = pass_box ? 9u : 0u;
    assert_ok_with_error(
        lammpc_session_nve_set_positions(session_ref, reference, 6, setter_box,
                                         setter_box_len),
        name);
    assert_ok_with_error(
        lammpc_session_nve_set_positions(session_img, image, 6, setter_box,
                                         setter_box_len),
        name);

    double pos_ref[6] = {0.0}, pos_img[6] = {0.0};
    LammpcNveObservation obs_ref = {0}, obs_img = {0};
    assert_ok_with_error(
        lammpc_session_nve_state(session_ref, pos_ref, NULL, 6, &obs_ref),
        name);
    assert_ok_with_error(
        lammpc_session_nve_state(session_img, pos_img, NULL, 6, &obs_img),
        name);
    assert_true(fabs(obs_ref.potential_energy_ev) > 1e-6);
    assert_true(
        fabs(obs_ref.potential_energy_ev - obs_img.potential_energy_ev) <
        1e-12);
    assert_positions_folded_and_matched(boundary, box9, pos_ref, pos_img, name);

    const double zero_vel[6] = {0.0};
    const double dt_ps = 1.0e-4;
    assert_ok_with_error(
        lammpc_session_nve_set_velocities(session_ref, zero_vel, 6), name);
    assert_ok_with_error(
        lammpc_session_nve_set_velocities(session_img, zero_vel, 6), name);
    assert_ok_with_error(lammpc_session_nve_step(session_ref, 1, dt_ps), name);
    assert_ok_with_error(lammpc_session_nve_step(session_img, 1, dt_ps), name);

    double pos1_ref[6] = {0.0}, pos1_img[6] = {0.0};
    double vel1_ref[6] = {0.0}, vel1_img[6] = {0.0};
    LammpcNveObservation obs1_ref = {0}, obs1_img = {0};
    assert_ok_with_error(
        lammpc_session_nve_state(session_ref, pos1_ref, vel1_ref, 6, &obs1_ref),
        name);
    assert_ok_with_error(
        lammpc_session_nve_state(session_img, pos1_img, vel1_img, 6, &obs1_img),
        name);

    const double force_scale =
        2.0 * kCuMassAmu * kMassVelocityToEv / (dt_ps * dt_ps);
    for (int i = 0; i < 6; ++i) {
      const double f_ref = force_scale * (pos1_ref[i] - pos_ref[i]);
      const double f_img = force_scale * (pos1_img[i] - pos_img[i]);
      assert_true(fabs(f_ref - f_img) < 1e-6);
      assert_true(fabs(vel1_ref[i] - vel1_img[i]) < 1e-12);
    }
    assert_true(
        fabs(obs1_ref.potential_energy_ev - obs1_img.potential_energy_ev) <
        1e-12);
    assert_true(
        fabs(obs1_ref.kinetic_energy_ev - obs1_img.kinetic_energy_ev) < 1e-12);

    lammpc_session_destroy(session_ref);
    lammpc_session_destroy(session_img);
  }
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

static void test_nve_set_positions_orthogonal_far_images(void **state) {
  (void)state;
  const double L = 7.992;
  const double box9[9] = {L, 0.0, 0.0, 0.0, L, 0.0, 0.0, 0.0, L};
  const double s0_ref[3] = {0.15, 0.25, 0.35};
  const double s1_ref[3] = {0.45, 0.55, 0.4};
  const double s0_img[3] = {4.15, -2.75, 5.35};
  const double s1_img[3] = {-4.55, 6.55, -3.6};
  double ref[6] = {0.0};
  double img[6] = {0.0};
  set_two_atom_positions(box9, s0_ref, s1_ref, ref);
  set_two_atom_positions(box9, s0_img, s1_img, img);
  run_lj_nve_set_positions_case("p p p", box9, ref, img,
                                "nve orthogonal far periodic image");

  const double x_image = nextafter(2.0 * L, 0.0);
  const double x_reference = x_image - L;
  const double ref_face[6] = {x_reference, 1.5, 0.75, 2.3, 2.5, 1.0};
  const double img_face[6] = {x_image, 1.5, 0.75, 2.3, 2.5, 1.0};
  run_lj_nve_set_positions_case("p p p", box9, ref_face, img_face,
                                "nve orthogonal nextafter(2L,0) image");
}

static void test_nve_set_positions_triclinic_far_images(void **state) {
  (void)state;
  const double box_many[9] = {12.0, 0.0, 0.0, 1.5, 11.5, 0.0, 0.75, 0.5, 11.0};
  const double s0_ref[3] = {0.125, 0.25, 0.5};
  const double s1_ref[3] = {0.35, 0.4, 0.6};
  const double s0_img[3] = {10.125, -8.75, 5.5};
  const double s1_img[3] = {-7.65, 12.4, -9.4};
  double ref[6] = {0.0};
  double img[6] = {0.0};
  set_two_atom_positions(box_many, s0_ref, s1_ref, ref);
  set_two_atom_positions(box_many, s0_img, s1_img, img);
  run_lj_nve_set_positions_case("p p p", box_many, ref, img,
                                "nve triclinic many-period image");

  const double box_plus2a[9] = {8.0, 0.0, 0.0, 1.0, 8.0, 0.0, 0.5, 0.25, 8.0};
  const double s1_plus2a[3] = {0.45, 0.2, 0.35};
  const double s0_plus2a_ref[3] = {0.0, 0.0, 0.3};
  const double s0_plus2a_img[3] = {2.0, 0.0, 0.3};
  set_two_atom_positions(box_plus2a, s0_plus2a_ref, s1_plus2a, ref);
  set_two_atom_positions(box_plus2a, s0_plus2a_img, s1_plus2a, img);
  run_lj_nve_set_positions_case("p p p", box_plus2a, ref, img,
                                "nve triclinic +2a face image");

  const double box_minus1a[9] = {
      10.86, 0.0, 0.0, 0.6, 10.86, 0.0, 0.3, 0.2, 10.86,
  };
  const double s1_minus1a[3] = {0.25, 0.35, 0.4};
  const double s0_minus1a_ref[3] = {0.0, 1.0, 0.1};
  const double s0_minus1a_img[3] = {-1.0, 1.0, 0.1};
  set_two_atom_positions(box_minus1a, s0_minus1a_ref, s1_minus1a, ref);
  set_two_atom_positions(box_minus1a, s0_minus1a_img, s1_minus1a, img);
  run_lj_nve_set_positions_case("p p p", box_minus1a, ref, img,
                                "nve triclinic -1a face image");
}

static void test_nve_set_positions_partly_periodic_far_images(void **state) {
  (void)state;
  const double box_tri[9] = {12.0, 0.0, 0.0, 1.5, 11.5, 0.0, 0.75, 0.5, 11.0};
  const double box_ortho[9] = {12.0, 0.0, 0.0, 0.0, 11.5, 0.0, 0.0, 0.0, 11.0};
  const double *boxes[2] = {box_ortho, box_tri};

  struct partial_case {
    const char *boundary;
    int periodic[3];
    int nonperiodic_axis;
    double nonperiodic_s0;
    double nonperiodic_s1;
  };
  const struct partial_case cases[] = {
      {"p p f", {1, 1, 0}, 2, 0.35, 0.55},
      {"p f p", {1, 0, 1}, 1, 0.35, 0.55},
      {"f p p", {0, 1, 1}, 0, 0.35, 0.45},
      {"p p s", {1, 1, 0}, 2, 1.15, 0.85},
      {"p s p", {1, 0, 1}, 1, -0.15, 0.15},
      {"s p p", {0, 1, 1}, 0, 1.15, 0.85},
  };

  for (size_t b = 0; b < 2; ++b) {
    const double *box9 = boxes[b];
    for (size_t c = 0; c < sizeof(cases) / sizeof(cases[0]); ++c) {
      const struct partial_case *item = &cases[c];
      double s0_ref[3] = {0.2, 0.3, 0.4};
      double s1_ref[3] = {0.45, 0.5, 0.25};
      s0_ref[item->nonperiodic_axis] = item->nonperiodic_s0;
      s1_ref[item->nonperiodic_axis] = item->nonperiodic_s1;
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
      run_lj_nve_set_positions_case(item->boundary, box9, ref, img,
                                    item->boundary);
    }
  }
}

int main(void) {
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_nve_energy_and_time_reversibility),
      cmocka_unit_test(test_nve_set_positions_orthogonal_far_images),
      cmocka_unit_test(test_nve_set_positions_triclinic_far_images),
      cmocka_unit_test(test_nve_set_positions_partly_periodic_far_images),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
