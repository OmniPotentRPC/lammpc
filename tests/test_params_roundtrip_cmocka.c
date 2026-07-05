/**
 * Cap'n Proto plumbing round-trips without any engine: ForceInput geometry
 * decode (with unit conversion), LammpsParams field decode, and the
 * PotentialResult writer.
 */
#include "lammpc_params.h"

#include <setjmp.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include <cmocka.h>

static size_t write_flat(struct capn *arena, capn_ptr root_struct,
                         unsigned char *dst, size_t dst_size) {
  capn_ptr root = capn_root(arena);
  assert_int_equal(capn_setp(root, 0, root_struct), 0);
  int written = capn_write_mem(arena, dst, dst_size, 0);
  assert_true(written > 0);
  return (size_t)written;
}

static void test_force_input_geometry_round_trip(void **state) {
  (void)state;
  struct capn arena;
  capn_init_malloc(&arena);
  ForceInput_ptr input = new_ForceInput(capn_root(&arena).seg);
  struct capn_segment *seg = capn_root(&arena).seg;

  capn_list64 pos = capn_new_list64(seg, 6);
  const double coords_bohr[6] = {0.0, 0.0, 0.0, 2.0, 0.0, 0.0};
  for (int i = 0; i < 6; ++i)
    capn_set64(pos, i, capn_from_f64(coords_bohr[i]));
  capn_list32 atmnrs = capn_new_list32(seg, 2);
  capn_set32(atmnrs, 0, 29);
  capn_set32(atmnrs, 1, 29);
  capn_list64 box = capn_new_list64(seg, 9);
  for (int i = 0; i < 9; ++i)
    capn_set64(box, i, capn_from_f64(i % 4 == 0 ? 20.0 : 0.0));

  struct ForceInput view;
  memset(&view, 0, sizeof(view));
  view.pos = pos;
  view.atmnrs = atmnrs;
  view.box = box;
  view.lengthUnit = (capn_text){4, "bohr", NULL};
  view.energyUnit = (capn_text){2, "eV", NULL};
  write_ForceInput(&view, input);

  unsigned char buf[4096];
  size_t written = write_flat(&arena, input.p, buf, sizeof(buf));
  capn_free(&arena);

  struct capn parse_arena;
  ForceInput_ptr parsed;
  assert_int_equal(lammpc_force_input_root(buf, written, &parse_arena,
                                           &parsed),
                   0);
  size_t n_atoms = 0;
  int has_cell = 0;
  assert_int_equal(lammpc_force_input_atom_count(parsed, &n_atoms, &has_cell),
                   0);
  assert_int_equal((int)n_atoms, 2);
  assert_int_equal(has_cell, 1);

  double positions_ang[6];
  int atomic_numbers[2];
  double cell_ang[9];
  assert_int_equal(lammpc_force_input_copy_geometry(parsed, positions_ang,
                                                    atomic_numbers, 2,
                                                    cell_ang, &has_cell),
                   0);
  assert_int_equal(atomic_numbers[0], 29);
  /* 2 bohr -> Angstrom. */
  assert_true(positions_ang[3] > 1.0583 && positions_ang[3] < 1.0584);
  assert_true(cell_ang[0] > 10.583 && cell_ang[0] < 10.584);

  double energy_factor = 0.0;
  double force_factor = 0.0;
  assert_int_equal(lammpc_force_input_result_factors(parsed, &energy_factor,
                                                     &force_factor),
                   0);
  /* Hartree -> eV. */
  assert_true(energy_factor > 27.211 && energy_factor < 27.212);
  lammpc_params_release(&parse_arena);
}

static void test_lammps_params_decode(void **state) {
  (void)state;
  struct capn arena;
  capn_init_malloc(&arena);
  struct capn_segment *seg = capn_root(&arena).seg;
  LammpsParams_ptr params = new_LammpsParams(seg);

  capn_ptr coeffs = capn_new_ptr_list(seg, 1);
  capn_text coeff = {17, "* * Cu_u3.eam.pot", NULL};
  assert_int_equal(capn_set_text(coeffs, 0, coeff), 0);
  capn_list32 type_map = capn_new_list32(seg, 1);
  capn_set32(type_map, 0, 29);

  struct LammpsParams view;
  memset(&view, 0, sizeof(view));
  view.unitsStyle = (capn_text){5, "metal", NULL};
  view.atomStyle = (capn_text){6, "atomic", NULL};
  view.pairStyle = (capn_text){9, "eam/alloy", NULL};
  view.pairCoeffs = coeffs;
  view.typeToAtomicNumber = type_map;
  view.newtonPair = 1;
  view.boundary = (capn_text){5, "p p p", NULL};
  write_LammpsParams(&view, params);

  unsigned char buf[4096];
  size_t written = write_flat(&arena, params.p, buf, sizeof(buf));
  capn_free(&arena);

  struct capn parse_arena;
  LammpsParams_ptr parsed;
  assert_int_equal(lammpc_params_root(buf, written, &parse_arena, &parsed), 0);
  struct LammpsParams decoded;
  read_LammpsParams(&decoded, parsed);
  assert_int_equal(decoded.pairStyle.len, 9);
  assert_memory_equal(decoded.pairStyle.str, "eam/alloy", 9);
  capn_resolve(&decoded.typeToAtomicNumber.p);
  assert_int_equal(decoded.typeToAtomicNumber.p.len, 1);
  assert_int_equal((int32_t)capn_get32(decoded.typeToAtomicNumber, 0), 29);
  capn_resolve(&decoded.pairCoeffs);
  assert_int_equal(decoded.pairCoeffs.len, 1);
  lammpc_params_release(&parse_arena);
}

static void test_potential_result_write_and_read(void **state) {
  (void)state;
  const double forces[6] = {0.1, 0.2, 0.3, -0.1, -0.2, -0.3};
  size_t need = lammpc_potential_result_flat_size(6);
  assert_true(need > 0);
  unsigned char *buf = (unsigned char *)malloc(need);
  assert_non_null(buf);
  size_t wrote = 0;
  assert_int_equal(lammpc_potential_result_write(-1.5, forces, 6, buf, need,
                                                 &wrote),
                   0);
  assert_true(wrote > 0);

  struct capn arena;
  assert_int_equal(capn_init_mem(&arena, buf, wrote, 0), 0);
  PotentialResult_ptr result;
  result.p = capn_getp(capn_root(&arena), 0, 1);
  assert_int_equal(result.p.type, CAPN_STRUCT);
  struct PotentialResult view;
  read_PotentialResult(&view, result);
  assert_true(view.energy == -1.5);
  capn_resolve(&view.forces.p);
  assert_int_equal(view.forces.p.len, 6);
  assert_true(capn_to_f64(capn_get64(view.forces, 2)) == 0.3);
  capn_free(&arena);
  free(buf);
}

int main(void) {
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_force_input_geometry_round_trip),
      cmocka_unit_test(test_lammps_params_decode),
      cmocka_unit_test(test_potential_result_write_and_read),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
