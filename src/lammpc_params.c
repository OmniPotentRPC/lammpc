#include "lammpc_params.h"

#include <limits.h>
#include <stdint.h>
#include <string.h>

static int text_eq_ci(capn_text text, const char *expected) {
  size_t expected_len = strlen(expected);
  if (text.len != (int)expected_len || !text.str)
    return 0;
  for (size_t i = 0; i < expected_len; ++i) {
    char a = text.str[i];
    char b = expected[i];
    if (a >= 'A' && a <= 'Z')
      a = (char)(a + ('a' - 'A'));
    if (b >= 'A' && b <= 'Z')
      b = (char)(b + ('a' - 'A'));
    if (a != b)
      return 0;
  }
  return 1;
}

static int list64_len(capn_list64 *list) {
  capn_resolve(&list->p);
  if (list->p.type == CAPN_NULL)
    return 0;
  if (list->p.type != CAPN_LIST)
    return -1;
  return list->p.len;
}

static int list32_len(capn_list32 *list) {
  capn_resolve(&list->p);
  if (list->p.type == CAPN_NULL)
    return 0;
  if (list->p.type != CAPN_LIST)
    return -1;
  return list->p.len;
}

static int force_input_length_factor(capn_text unit, double *factor) {
  if (!factor)
    return -1;
  if (!unit.str || unit.len <= 0 || text_eq_ci(unit, "angstrom") ||
      text_eq_ci(unit, "angstroms") || text_eq_ci(unit, "a") ||
      text_eq_ci(unit, "aa")) {
    *factor = 1.0;
    return 0;
  }
  if (text_eq_ci(unit, "bohr") || text_eq_ci(unit, "au") ||
      text_eq_ci(unit, "a.u.") || text_eq_ci(unit, "atomic")) {
    *factor = LAMMPC_BOHR_TO_ANGSTROM;
    return 0;
  }
  if (text_eq_ci(unit, "nm") || text_eq_ci(unit, "nanometer")) {
    *factor = 10.0;
    return 0;
  }
  return -1;
}

static int force_input_energy_factor(capn_text unit, double *factor) {
  if (!factor)
    return -1;
  if (!unit.str || unit.len <= 0 || text_eq_ci(unit, "hartree") ||
      text_eq_ci(unit, "ha") || text_eq_ci(unit, "au") ||
      text_eq_ci(unit, "a.u.")) {
    *factor = 1.0;
    return 0;
  }
  if (text_eq_ci(unit, "ev") || text_eq_ci(unit, "electronvolt")) {
    *factor = LAMMPC_HARTREE_TO_EV;
    return 0;
  }
  if (text_eq_ci(unit, "ry") || text_eq_ci(unit, "rydberg")) {
    *factor = 2.0;
    return 0;
  }
  if (text_eq_ci(unit, "kj/mol") || text_eq_ci(unit, "kjmol")) {
    *factor = 2625.499639479;
    return 0;
  }
  if (text_eq_ci(unit, "kcal/mol") || text_eq_ci(unit, "kcalmol")) {
    *factor = 627.5094740631;
    return 0;
  }
  return -1;
}

static int flat_root_struct(const void *msg, size_t msg_size_bytes,
                            struct capn *arena, capn_ptr *root_struct) {
  if (!msg || msg_size_bytes == 0 || !arena || !root_struct)
    return -1;
  memset(arena, 0, sizeof(*arena));
  if (capn_init_mem(arena, (const uint8_t *)msg, msg_size_bytes, 0) != 0)
    return -1;
  *root_struct = capn_getp(capn_root(arena), 0, 1);
  if (root_struct->type != CAPN_STRUCT) {
    lammpc_params_release(arena);
    return -1;
  }
  return 0;
}

int lammpc_force_input_root(const void *force_input_capnp,
                            size_t force_input_capnp_size_bytes,
                            struct capn *arena, ForceInput_ptr *force_input) {
  if (!force_input)
    return -1;
  memset(force_input, 0, sizeof(*force_input));
  return flat_root_struct(force_input_capnp, force_input_capnp_size_bytes,
                          arena, &force_input->p);
}

int lammpc_params_root(const void *params_capnp,
                       size_t params_capnp_size_bytes, struct capn *arena,
                       LammpsParams_ptr *params) {
  if (!params)
    return -1;
  memset(params, 0, sizeof(*params));
  return flat_root_struct(params_capnp, params_capnp_size_bytes, arena,
                          &params->p);
}

void lammpc_params_release(struct capn *arena) {
  if (!arena)
    return;
  capn_free(arena);
  memset(arena, 0, sizeof(*arena));
}

int lammpc_force_input_atom_count(ForceInput_ptr force_input, size_t *n_atoms,
                                  int *has_cell) {
  if (force_input.p.type == CAPN_NULL || !n_atoms || !has_cell)
    return -1;
  struct ForceInput view;
  read_ForceInput(&view, force_input);
  int n_pos = list64_len(&view.pos);
  int n_z = list32_len(&view.atmnrs);
  int n_box = list64_len(&view.box);
  if (n_pos < 0 || n_z <= 0 || n_box < 0)
    return -1;
  if (n_pos != n_z * 3)
    return -1;
  if (n_box != 0 && n_box != 9)
    return -1;
  *n_atoms = (size_t)n_z;
  *has_cell = n_box == 9 ? 1 : 0;
  return 0;
}

int lammpc_force_input_copy_geometry(ForceInput_ptr force_input,
                                     double *positions_ang,
                                     int *atomic_numbers, size_t atom_capacity,
                                     double *cell_ang, int *has_cell) {
  if (force_input.p.type == CAPN_NULL || !positions_ang || !atomic_numbers ||
      !cell_ang || !has_cell)
    return -1;
  struct ForceInput view;
  read_ForceInput(&view, force_input);
  size_t n_atoms = 0;
  int local_has_cell = 0;
  if (lammpc_force_input_atom_count(force_input, &n_atoms, &local_has_cell) !=
      0)
    return -1;
  if (atom_capacity < n_atoms)
    return -1;
  double length_factor = 1.0;
  if (force_input_length_factor(view.lengthUnit, &length_factor) != 0)
    return -1;
  capn_resolve(&view.pos.p);
  capn_resolve(&view.atmnrs.p);
  capn_resolve(&view.box.p);
  for (size_t i = 0; i < n_atoms; ++i)
    atomic_numbers[i] = (int)(int32_t)capn_get32(view.atmnrs, (int)i);
  for (size_t i = 0; i < n_atoms * 3u; ++i)
    positions_ang[i] =
        capn_to_f64(capn_get64(view.pos, (int)i)) * length_factor;
  for (size_t i = 0; i < 9; ++i)
    cell_ang[i] =
        local_has_cell
            ? capn_to_f64(capn_get64(view.box, (int)i)) * length_factor
            : 0.0;
  *has_cell = local_has_cell;
  return 0;
}

int lammpc_force_input_result_factors(ForceInput_ptr force_input,
                                      double *energy_factor,
                                      double *force_factor) {
  if (force_input.p.type == CAPN_NULL || !energy_factor || !force_factor)
    return -1;
  struct ForceInput view;
  read_ForceInput(&view, force_input);
  double length_factor = 1.0;
  double energy = 1.0;
  if (force_input_length_factor(view.lengthUnit, &length_factor) != 0 ||
      force_input_energy_factor(view.energyUnit, &energy) != 0)
    return -1;
  *energy_factor = energy;
  *force_factor = energy * length_factor / LAMMPC_BOHR_TO_ANGSTROM;
  return 0;
}

size_t lammpc_potential_result_flat_size(size_t force_count) {
  /* One natoms*3 Float64 forces list plus capnp framing headroom. */
  const size_t overhead = 512u;
  if (force_count > (SIZE_MAX - overhead) / 8u)
    return 0;
  return overhead + force_count * 8u;
}

int lammpc_potential_result_write(double energy, const double *forces,
                                  size_t force_count,
                                  void *potential_result_capnp,
                                  size_t potential_result_capacity_bytes,
                                  size_t *potential_result_size_bytes) {
  if (!forces || !potential_result_capnp || !potential_result_size_bytes ||
      force_count > (size_t)INT_MAX)
    return -1;
  size_t required = lammpc_potential_result_flat_size(force_count);
  *potential_result_size_bytes = required;
  if (required == 0 || potential_result_capacity_bytes < required)
    return -1;

  struct capn arena;
  capn_init_malloc(&arena);
  capn_ptr root = capn_root(&arena);
  if (root.type == CAPN_NULL) {
    capn_free(&arena);
    return -1;
  }
  PotentialResult_ptr result = new_PotentialResult(root.seg);
  capn_list64 force_list = capn_new_list64(root.seg, (int)force_count);
  if (result.p.type == CAPN_NULL ||
      (force_count > 0 && force_list.p.type == CAPN_NULL)) {
    capn_free(&arena);
    return -1;
  }
  for (size_t i = 0; i < force_count; ++i)
    capn_set64(force_list, (int)i, capn_from_f64(forces[i]));

  struct PotentialResult view;
  memset(&view, 0, sizeof(view));
  view.energy = energy;
  view.forces = force_list;
  view.embedMdPropsSkipped = 1;
  write_PotentialResult(&view, result);
  if (capn_setp(root, 0, result.p) != 0) {
    capn_free(&arena);
    return -1;
  }
  int written = capn_write_mem(&arena, (uint8_t *)potential_result_capnp,
                               potential_result_capacity_bytes, 0);
  capn_free(&arena);
  if (written < 0)
    return -1;
  *potential_result_size_bytes = (size_t)written;
  return 0;
}
