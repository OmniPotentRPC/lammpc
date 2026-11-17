#define _GNU_SOURCE

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

static size_t make_config(const char *pair_style, const char *pair_coeff_line,
                          int atomic_number, double mass,
                          const char *boundary, unsigned char *dst,
                          size_t dst_size) {
  struct capn arena;
  capn_init_malloc(&arena);
  struct capn_segment *seg = capn_root(&arena).seg;

  LammpsParams_ptr params = new_LammpsParams(seg);
  capn_ptr coeffs = capn_new_ptr_list(seg, 1);
  capn_text coeff = {(int)strlen(pair_coeff_line), pair_coeff_line, 0};
  assert_int_equal(capn_set_text(coeffs, 0, coeff), 0);
  capn_list32 type_map = capn_new_list32(seg, 1);
  capn_set32(type_map, 0, atomic_number);
  capn_list64 masses = capn_new_list64(seg, 1);
  capn_set64(masses, 0, capn_from_f64(mass));

  struct LammpsParams pview;
  memset(&pview, 0, sizeof(pview));
  pview.unitsStyle = (capn_text){5, "metal", NULL};
  pview.atomStyle = (capn_text){6, "atomic", NULL};
  pview.pairStyle = (capn_text){(int)strlen(pair_style), pair_style, NULL};
  pview.pairCoeffs = coeffs;
  pview.typeToAtomicNumber = type_map;
  pview.masses = masses;
  pview.newtonPair = 1;
  pview.boundary =
      (capn_text){(int)strlen(boundary), (char *)boundary, NULL};
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
                        const double *box9, int with_box, unsigned char *dst,
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
  if (with_box) {
    capn_list64 box = capn_new_list64(seg, 9);
    for (int i = 0; i < 9; ++i)
      capn_set64(box, i, capn_from_f64(box9[i]));
    view.box = box;
  }
  view.lengthUnit = (capn_text){8, "angstrom", NULL};
  view.energyUnit = (capn_text){2, "eV", NULL};
  write_ForceInput(&view, input);

  size_t written = write_flat(&arena, input.p, dst, dst_size);
  capn_free(&arena);
  return written;
}

static void decode_result(const unsigned char *buf, size_t size, double *energy,
                          double *forces, int n_atoms, double stress9[9]) {
  struct capn arena;
  assert_int_equal(capn_init_mem(&arena, buf, size, 0), 0);
  PotentialResult_ptr result;
  result.p = capn_getp(capn_root(&arena), 0, 1);
  assert_int_equal(result.p.type, CAPN_STRUCT);
  struct PotentialResult view;
  read_PotentialResult(&view, result);
  *energy = view.energy;
  capn_resolve(&view.forces.p);
  assert_int_equal(view.forces.p.len, n_atoms * 3);
  for (int i = 0; i < n_atoms * 3; ++i)
    forces[i] = capn_to_f64(capn_get64(view.forces, i));
  capn_resolve(&view.stress.p);
  assert_int_equal(view.stress.p.len, 9);
  for (int i = 0; i < 9; ++i)
    stress9[i] = capn_to_f64(capn_get64(view.stress, i));
  capn_free(&arena);
}

static int parse_tagged_line(const char *line, const char *tag, double *values,
                             int n_values) {
  size_t tag_len = strlen(tag);
  if (strncmp(line, tag, tag_len) != 0)
    return 0;
  const char *p = line + tag_len;
  for (int i = 0; i < n_values; ++i) {
    char *end = NULL;
    values[i] = strtod(p, &end);
    if (end == p)
      return -1;
    p = end;
  }
  return 1;
}

static void run_lmp_reference(const char *lammps_root, const char *lmp_exe,
                              const char *pair_style,
                              const char *pair_coeff, const char *mass_line,
                              const double *positions_ang, int n_atoms,
                              const double *box9, double *energy,
                              double *forces, double stress6[6]) {
  char tmp_template[] = "lammpc-ref-XXXXXX";
  char *tmp_dir = mkdtemp(tmp_template);
  assert_non_null(tmp_dir);

  char input_path[1024];
  char dump_path[1024];
  const char *dump_file = "dump.forces";
  snprintf(input_path, sizeof(input_path), "%s/in.lmp", tmp_dir);
  snprintf(dump_path, sizeof(dump_path), "%s/%s", tmp_dir, dump_file);

  FILE *in = fopen(input_path, "w");
  assert_non_null(in);
  fprintf(in, "units metal\n");
  fprintf(in, "atom_style atomic\n");
  fprintf(in, "boundary p p p\n");
  fprintf(in, "atom_modify map array\n");
  fprintf(in, "region box prism 0 %.16g 0 %.16g 0 %.16g %.16g %.16g %.16g "
              "units box\n",
          box9[0], box9[4], box9[8], box9[3], box9[6], box9[7]);
  fprintf(in, "create_box 1 box\n");
  for (int i = 0; i < n_atoms; ++i) {
    fprintf(in, "create_atoms 1 single %.16g %.16g %.16g units box\n",
            positions_ang[3 * i + 0], positions_ang[3 * i + 1],
            positions_ang[3 * i + 2]);
  }
  fprintf(in, "mass 1 %s\n", mass_line);
  fprintf(in, "pair_style %s\n", pair_style);
  fprintf(in, "pair_coeff %s\n", pair_coeff);
  fprintf(in, "compute p all pressure NULL virial\n");
  fprintf(in, "variable pe equal pe\n");
  fprintf(in, "variable sxx equal -c_p[1]/1602176.6208\n");
  fprintf(in, "variable syy equal -c_p[2]/1602176.6208\n");
  fprintf(in, "variable szz equal -c_p[3]/1602176.6208\n");
  fprintf(in, "variable sxy equal -c_p[4]/1602176.6208\n");
  fprintf(in, "variable sxz equal -c_p[5]/1602176.6208\n");
  fprintf(in, "variable syz equal -c_p[6]/1602176.6208\n");
  fprintf(in, "dump d all custom 1 %s id fx fy fz\n", dump_file);
  fprintf(in, "dump_modify d sort id\n");
  fprintf(in, "run 0\n");
  fprintf(in, "print \"__LAMMPS_PE ${pe}\"\n");
  fprintf(in, "print \"__LAMMPS_STRESS ${sxx} ${syy} ${szz} ${sxy} ${sxz} ${syz}\"\n");
  fclose(in);

  char cmd[2048];
  snprintf(cmd, sizeof(cmd),
           "cd '%s' && LD_LIBRARY_PATH='%s/lib' '%s' -log none -in in.lmp 2>&1",
           tmp_dir, lammps_root, lmp_exe);
  FILE *pipe = popen(cmd, "r");
  assert_non_null(pipe);
  char line[2048];
  char captured[8192];
  captured[0] = '\0';
  int got_pe = 0;
  int got_stress = 0;
  while (fgets(line, sizeof(line), pipe) != NULL) {
    if (strlen(captured) + strlen(line) + 1 < sizeof(captured))
      strcat(captured, line);
    double vals[6] = {0};
    int rc = parse_tagged_line(line, "__LAMMPS_PE ", vals, 1);
    if (rc == 1) {
      *energy = vals[0];
      got_pe = 1;
      continue;
    }
    rc = parse_tagged_line(line, "__LAMMPS_STRESS ", vals, 6);
    if (rc == 1) {
      for (int i = 0; i < 6; ++i)
        stress6[i] = vals[i];
      got_stress = 1;
    }
  }
  int status = pclose(pipe);
  if (status != 0)
    fail_msg("lmp failed with status=%d output=\n%s", status, captured);
  assert_true(got_pe);
  assert_true(got_stress);

  FILE *dump = fopen(dump_path, "r");
  assert_non_null(dump);
  int found_atoms = 0;
  int atoms_to_read = 0;
  while (fgets(line, sizeof(line), dump) != NULL) {
    if (strncmp(line, "ITEM: NUMBER OF ATOMS", 21) == 0) {
      assert_non_null(fgets(line, sizeof(line), dump));
      atoms_to_read = atoi(line);
      assert_int_equal(atoms_to_read, n_atoms);
    } else if (strncmp(line, "ITEM: ATOMS", 11) == 0) {
      for (int i = 0; i < atoms_to_read; ++i) {
        int id = 0;
        double fx = 0, fy = 0, fz = 0;
        assert_non_null(fgets(line, sizeof(line), dump));
        assert_int_equal(sscanf(line, "%d %lf %lf %lf", &id, &fx, &fy, &fz), 4);
        assert_true(id >= 1 && id <= n_atoms);
        forces[3 * (id - 1) + 0] = fx;
        forces[3 * (id - 1) + 1] = fy;
        forces[3 * (id - 1) + 2] = fz;
      }
      found_atoms = 1;
    }
  }
  fclose(dump);
  assert_true(found_atoms);
}

static void eval_lammpc(const unsigned char *config, size_t config_size,
                        const unsigned char *step, size_t step_size,
                        int n_atoms, double *energy, double *forces,
                        double stress9[9]) {
  LammpcSession *session =
      lammpc_session_create_from_config(config, config_size);
  assert_non_null(session);
  size_t need = lammpc_potential_result_size_for_force_input(step, step_size);
  assert_true(need > 0);
  unsigned char *out = (unsigned char *)malloc(need);
  assert_non_null(out);
  size_t wrote = 0;
  LammpcResult r = lammpc_session_calculate_result(session, step, step_size,
                                                   out, need, &wrote);
  if (!r.ok)
    fail_msg("lammpc eval failed: %s", r.message);
  decode_result(out, wrote, energy, forces, n_atoms, stress9);
  free(out);
  lammpc_session_destroy(session);
}

static void eval_lammpc_raw(const unsigned char *config, size_t config_size,
                            const double *positions, const int *z, int n_atoms,
                            const double box9[9], double *energy,
                            double *forces, double stress9[9]) {
  const double ev_per_hartree = 27.211386245988;
  const double angstrom_per_bohr = 0.529177210903;
  LammpcSession *session =
      lammpc_session_create_from_config(config, config_size);
  assert_non_null(session);
  double *gradient = calloc((size_t)n_atoms * 3, sizeof(double));
  assert_non_null(gradient);
  double raw_stress[9] = {0};
  LammpcResult result = lammpc_session_energy_gradient(
      session, n_atoms, positions, z, box9, 9, gradient, raw_stress);
  if (!result.ok)
    fail_msg("lammpc raw eval failed: %s", result.message);
  *energy = result.energy_h * ev_per_hartree;
  for (int i = 0; i < n_atoms * 3; ++i)
    forces[i] = -gradient[i] * ev_per_hartree / angstrom_per_bohr;
  const double stress_factor =
      ev_per_hartree / (angstrom_per_bohr * angstrom_per_bohr *
                        angstrom_per_bohr);
  for (int i = 0; i < 9; ++i)
    stress9[i] = raw_stress[i] * stress_factor;
  free(gradient);
  lammpc_session_destroy(session);
}

static void finite_difference_force_check(const unsigned char *config,
                                          size_t config_size, double *positions,
                                          const int *z, int n_atoms,
                                          const double *box9) {
  const double h = 1e-4;
  unsigned char step0[8192];
  size_t step0_size =
      make_step(positions, z, n_atoms, box9, 1, step0, sizeof(step0));
  double e0 = 0.0, f0[24] = {0}, s0[9] = {0};
  eval_lammpc(config, config_size, step0, step0_size, n_atoms, &e0, f0, s0);

  positions[0] += h;
  unsigned char step_p[8192];
  size_t step_p_size =
      make_step(positions, z, n_atoms, box9, 1, step_p, sizeof(step_p));
  double ep = 0.0, fp[24] = {0}, sp[9] = {0};
  eval_lammpc(config, config_size, step_p, step_p_size, n_atoms, &ep, fp, sp);

  positions[0] -= 2.0 * h;
  unsigned char step_m[8192];
  size_t step_m_size =
      make_step(positions, z, n_atoms, box9, 1, step_m, sizeof(step_m));
  double em = 0.0, fm[24] = {0}, sm[9] = {0};
  eval_lammpc(config, config_size, step_m, step_m_size, n_atoms, &em, fm, sm);
  positions[0] += h;

  const double fd_force = -(ep - em) / (2.0 * h);
  assert_true(fabs(fd_force - f0[0]) < 2e-3);
}

static void check_case(const char *pair_style, const char *pair_coeff,
                       int z, double mass, const double *positions, int n_atoms,
                       const double *box9, int run_fd_check) {
  const char *lammps_root = getenv("LAMMPS_ROOT");
  assert_non_null(lammps_root);
  char lmp_exe[1024];
  snprintf(lmp_exe, sizeof(lmp_exe), "%s/bin/lmp", lammps_root);
  assert_int_equal(access(lmp_exe, X_OK), 0);

  unsigned char config[16384];
  size_t config_size = make_config(pair_style, pair_coeff, z, mass, "p p p",
                                   config, sizeof(config));
  int zs[8];
  for (int i = 0; i < n_atoms; ++i)
    zs[i] = z;
  unsigned char step[8192];
  size_t step_size =
      make_step(positions, zs, n_atoms, box9, 1, step, sizeof(step));

  double e_lammpc = 0.0, f_lammpc[24] = {0}, s_lammpc[9] = {0};
  eval_lammpc(config, config_size, step, step_size, n_atoms, &e_lammpc,
              f_lammpc, s_lammpc);
  double e_raw = 0.0, f_raw[24] = {0}, s_raw[9] = {0};
  eval_lammpc_raw(config, config_size, positions, zs, n_atoms, box9, &e_raw,
                  f_raw, s_raw);
  assert_true(fabs(e_raw - e_lammpc) < 1e-10);
  for (int i = 0; i < 3 * n_atoms; ++i)
    assert_true(fabs(f_raw[i] - f_lammpc[i]) < 1e-10);
  for (int i = 0; i < 9; ++i)
    assert_true(fabs(s_raw[i] - s_lammpc[i]) < 1e-10);

  double e_lmp = 0.0, f_lmp[24] = {0}, s_lmp[6] = {0};
  char mass_line[64];
  snprintf(mass_line, sizeof(mass_line), "%.16g", mass);
  run_lmp_reference(lammps_root, lmp_exe, pair_style, pair_coeff, mass_line,
                    positions, n_atoms, box9, &e_lmp, f_lmp, s_lmp);

  assert_true(fabs(e_lammpc - e_lmp) < 1e-8);
  double max_force_delta = 0.0;
  for (int i = 0; i < 3 * n_atoms; ++i) {
    double delta = fabs(f_lammpc[i] - f_lmp[i]);
    if (delta > max_force_delta)
      max_force_delta = delta;
  }
  assert_true(max_force_delta < 2e-6);
  assert_true(fabs(s_lammpc[0] - s_lmp[0]) < 3e-8);
  assert_true(fabs(s_lammpc[4] - s_lmp[1]) < 3e-8);
  assert_true(fabs(s_lammpc[8] - s_lmp[2]) < 3e-8);
  assert_true(fabs(s_lammpc[1] - s_lmp[3]) < 3e-8);
  assert_true(fabs(s_lammpc[2] - s_lmp[4]) < 3e-8);
  assert_true(fabs(s_lammpc[5] - s_lmp[5]) < 3e-8);

  if (run_fd_check) {
    double pos_copy[24];
    memcpy(pos_copy, positions, sizeof(double) * 3 * n_atoms);
    finite_difference_force_check(config, config_size, pos_copy, zs, n_atoms,
                                  box9);
  }
}

static void test_eam_alloy_cu_vs_lmp(void **state) {
  (void)state;
  const char *lammps_root = getenv("LAMMPS_ROOT");
  assert_non_null(lammps_root);
  char potential[1024];
  snprintf(potential, sizeof(potential),
           "%s/share/lammps/potentials/Cu_mishin1.eam.alloy",
           lammps_root);
  assert_int_equal(access(potential, R_OK), 0);
  char coeff[1200];
  snprintf(coeff, sizeof(coeff), "* * %s Cu", potential);

  const double box9[9] = {7.22, 0.0, 0.0, 0.0, 7.22, 0.0, 0.0, 0.0, 7.22};
  const double pos[6] = {0.0, 0.0, 0.0, 1.805, 1.805, 1.805};
  check_case("eam/alloy", coeff, 29, 63.546, pos, 2, box9, 0);
}

static void test_tersoff_si_vs_lmp_triclinic(void **state) {
  (void)state;
  const char *lammps_root = getenv("LAMMPS_ROOT");
  assert_non_null(lammps_root);
  char potential[1024];
  snprintf(potential, sizeof(potential), "%s/share/lammps/potentials/Si.tersoff",
           lammps_root);
  assert_int_equal(access(potential, R_OK), 0);
  char coeff[1200];
  snprintf(coeff, sizeof(coeff), "* * %s Si", potential);

  const double box9[9] = {10.86, 0.0, 0.0, 0.6, 10.86, 0.0, 0.3, 0.2, 10.86};
  const double pos[6] = {0.0, 0.0, 0.0, 2.35, 2.35, 2.35};
  check_case("tersoff", coeff, 14, 28.0855, pos, 2, box9, 1);
}

int main(void) {
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_eam_alloy_cu_vs_lmp),
      cmocka_unit_test(test_tersoff_si_vs_lmp_triclinic),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
