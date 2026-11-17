#include "lammpc.h"
#include "lammpc_params.h"

#include "library.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

thread_local char g_last_error[512] = "";

void store_error(const std::string &message) {
  std::snprintf(g_last_error, sizeof(g_last_error), "%s", message.c_str());
}

void clear_error() { g_last_error[0] = '\0'; }

LammpcResult fail_result(const char *message) {
  LammpcResult r;
  r.ok = 0;
  r.energy_h = 0.0;
  std::snprintf(r.message, sizeof(r.message), "%s", message);
  return r;
}

std::string text_string(capn_text text) {
  return text.str && text.len > 0 ? std::string(text.str, (size_t)text.len)
                                  : std::string();
}

struct DecodedParams {
  std::string units_style = "metal";
  std::string atom_style = "atomic";
  std::string pair_style;
  std::vector<std::string> pair_coeffs;
  std::vector<int> type_to_z;
  std::vector<double> masses;
  std::string boundary = "p p p";
  std::string suffix;
  std::vector<std::string> extra_setup;
};

int decode_lammps_params(LammpsParams_ptr params, DecodedParams *out) {
  struct LammpsParams view;
  read_LammpsParams(&view, params);
  if (view.unitsStyle.len > 0)
    out->units_style = text_string(view.unitsStyle);
  if (view.atomStyle.len > 0)
    out->atom_style = text_string(view.atomStyle);
  out->pair_style = text_string(view.pairStyle);
  if (view.boundary.len > 0)
    out->boundary = text_string(view.boundary);
  out->suffix = text_string(view.suffix);

  capn_resolve(&view.pairCoeffs);
  if (view.pairCoeffs.type == CAPN_PTR_LIST ||
      (view.pairCoeffs.type == CAPN_LIST && view.pairCoeffs.len > 0)) {
    static const capn_text empty_text = {0, "", 0};
    for (int i = 0; i < view.pairCoeffs.len; ++i)
      out->pair_coeffs.push_back(
          text_string(capn_get_text(view.pairCoeffs, i, empty_text)));
  }
  capn_resolve(&view.typeToAtomicNumber.p);
  for (int i = 0; i < view.typeToAtomicNumber.p.len; ++i)
    out->type_to_z.push_back(
        (int)(int32_t)capn_get32(view.typeToAtomicNumber, i));
  capn_resolve(&view.masses.p);
  for (int i = 0; i < view.masses.p.len; ++i)
    out->masses.push_back(capn_to_f64(capn_get64(view.masses, i)));
  capn_resolve(&view.extraSetup);
  if (view.extraSetup.type == CAPN_PTR_LIST ||
      (view.extraSetup.type == CAPN_LIST && view.extraSetup.len > 0)) {
    static const capn_text empty_text = {0, "", 0};
    for (int i = 0; i < view.extraSetup.len; ++i)
      out->extra_setup.push_back(
          text_string(capn_get_text(view.extraSetup, i, empty_text)));
  }

  if (out->units_style != "metal") {
    store_error("only unitsStyle=metal is supported");
    return -1;
  }
  if (out->atom_style != "atomic") {
    store_error("only atomStyle=atomic is supported");
    return -1;
  }
  if (out->pair_style.empty()) {
    store_error("LammpsParams.pairStyle is required");
    return -1;
  }
  if (out->type_to_z.empty()) {
    store_error("LammpsParams.typeToAtomicNumber is required");
    return -1;
  }
  return 0;
}

constexpr double kEvPerHartree = 27.211386245988;
constexpr double kAngstromPerBohr = 0.529177210903;
constexpr double kBarPerEvA3 = 1602176.6208;

struct LammpcSession {
  void *lmp = nullptr;
  DecodedParams params;
  size_t n_atoms = 0;
  int has_cell = 0;
  double cell_ang[9] = {0};
  std::vector<int> atom_types;
};

int type_for_atomic_number(const DecodedParams &params, int z) {
  for (size_t t = 0; t < params.type_to_z.size(); ++t)
    if (params.type_to_z[t] == z)
      return (int)t + 1;
  return -1;
}

int run_cmd(void *lmp, const std::string &cmd) {
  (void)lammps_command(lmp, cmd.c_str());
  if (lammps_has_error(lmp)) {
    char buffer[512] = {0};
    if (lammps_get_last_error_message(lmp, buffer, (int)sizeof(buffer)) > 0)
      store_error(buffer);
    else
      store_error("LAMMPS command failed: " + cmd);
    return -1;
  }
  return 0;
}

void destroy_instance(LammpcSession *session) {
  if (!session || !session->lmp)
    return;
  lammps_close(session->lmp);
  session->lmp = nullptr;
}

int build_instance(LammpcSession *session, size_t n_atoms,
                   const double *positions_ang, const int *atomic_numbers,
                   const double *cell_ang, int has_cell) {
  destroy_instance(session);
  const DecodedParams &p = session->params;

  char arg0[] = "lammpc";
  char arg1[] = "-log";
  char arg2[] = "none";
  char arg3[] = "-screen";
  char arg4[] = "none";
  char *args[] = {arg0, arg1, arg2, arg3, arg4};
  void *lmp = lammps_open_no_mpi(5, args, nullptr);
  if (!lmp) {
    store_error("lammps_open_no_mpi failed");
    return -1;
  }

  if (run_cmd(lmp, "units metal") != 0 || run_cmd(lmp, "atom_style atomic") != 0 ||
      run_cmd(lmp, "atom_modify map array") != 0) {
    lammps_close(lmp);
    return -1;
  }
  if (run_cmd(lmp, "boundary " + p.boundary) != 0) {
    lammps_close(lmp);
    return -1;
  }

  std::string region_cmd;
  if (has_cell) {
    if (std::fabs(cell_ang[1]) > 1e-12 || std::fabs(cell_ang[2]) > 1e-12 ||
        std::fabs(cell_ang[5]) > 1e-12) {
      store_error("only restricted triclinic cells are supported");
      lammps_close(lmp);
      return -1;
    }
    if (!(cell_ang[0] > 0.0 && cell_ang[4] > 0.0 && cell_ang[8] > 0.0)) {
      store_error("cell diagonal entries must be positive");
      lammps_close(lmp);
      return -1;
    }
    region_cmd = "region box prism 0 " + std::to_string(cell_ang[0]) + " 0 " +
                 std::to_string(cell_ang[4]) + " 0 " +
                 std::to_string(cell_ang[8]) + " " + std::to_string(cell_ang[3]) +
                 " " + std::to_string(cell_ang[6]) + " " +
                 std::to_string(cell_ang[7]) + " units box";
  } else {
    double pad = 50.0;
    double mins[3] = {1e300, 1e300, 1e300};
    double maxs[3] = {-1e300, -1e300, -1e300};
    for (size_t i = 0; i < n_atoms; ++i) {
      for (int d = 0; d < 3; ++d) {
        const double v = positions_ang[3 * i + d];
        if (v < mins[d])
          mins[d] = v;
        if (v > maxs[d])
          maxs[d] = v;
      }
    }
    region_cmd = "region box block " + std::to_string(mins[0] - pad) + " " +
                 std::to_string(maxs[0] + pad) + " " +
                 std::to_string(mins[1] - pad) + " " +
                 std::to_string(maxs[1] + pad) + " " +
                 std::to_string(mins[2] - pad) + " " +
                 std::to_string(maxs[2] + pad) + " units box";
  }
  if (run_cmd(lmp, region_cmd) != 0) {
    lammps_close(lmp);
    return -1;
  }
  if (run_cmd(lmp, "create_box " + std::to_string((int)p.type_to_z.size()) +
                       " box") != 0) {
    lammps_close(lmp);
    return -1;
  }

  for (size_t t = 0; t < p.type_to_z.size(); ++t) {
    double mass = (t < p.masses.size()) ? p.masses[t] : 1.0;
    if (run_cmd(lmp, "mass " + std::to_string((int)t + 1) + " " +
                         std::to_string(mass)) != 0) {
      lammps_close(lmp);
      return -1;
    }
  }

  session->atom_types.assign(n_atoms, 0);
  for (size_t i = 0; i < n_atoms; ++i) {
    int type = type_for_atomic_number(p, atomic_numbers[i]);
    if (type < 1) {
      store_error("ForceInput atomic number missing from typeToAtomicNumber");
      lammps_close(lmp);
      return -1;
    }
    session->atom_types[i] = type;
    if (run_cmd(lmp, "create_atoms " + std::to_string(type) + " single " +
                         std::to_string(positions_ang[3 * i + 0]) + " " +
                         std::to_string(positions_ang[3 * i + 1]) + " " +
                         std::to_string(positions_ang[3 * i + 2]) +
                         " units box") != 0) {
      lammps_close(lmp);
      return -1;
    }
  }

  if (!p.suffix.empty() && run_cmd(lmp, "suffix " + p.suffix) != 0) {
    lammps_close(lmp);
    return -1;
  }
  if (run_cmd(lmp, "pair_style " + p.pair_style) != 0) {
    lammps_close(lmp);
    return -1;
  }
  for (const auto &coeff : p.pair_coeffs) {
    if (run_cmd(lmp, "pair_coeff " + coeff) != 0) {
      lammps_close(lmp);
      return -1;
    }
  }
  for (const auto &line : p.extra_setup) {
    if (!line.empty() && run_cmd(lmp, line) != 0) {
      lammps_close(lmp);
      return -1;
    }
  }
  if (run_cmd(lmp, "thermo_style custom step pe pxx pyy pzz pxy pxz pyz") != 0 ||
      run_cmd(lmp, "run 0") != 0) {
    lammps_close(lmp);
    return -1;
  }

  session->lmp = lmp;
  session->n_atoms = n_atoms;
  session->has_cell = has_cell;
  std::memcpy(session->cell_ang, cell_ang, sizeof(session->cell_ang));
  return 0;
}

bool geometry_matches(const LammpcSession *session, size_t n_atoms,
                      const double *cell_ang, int has_cell) {
  if (!session->lmp || session->n_atoms != n_atoms || session->has_cell != has_cell)
    return false;
  for (int i = 0; i < 9; ++i)
    if (session->cell_ang[i] != cell_ang[i])
      return false;
  return true;
}

int evaluate(LammpcSession *session, size_t n_atoms, const double *positions_ang,
             double *energy_h, double *forces_h_bohr, double *stress_h_bohr3) {
  if (!session || !session->lmp) {
    store_error("no live LAMMPS session");
    return -1;
  }
  lammps_scatter_atoms(session->lmp, "x", 1, 3, (void *)positions_ang);
  if (run_cmd(session->lmp, "run 0 pre yes post no") != 0)
    return -1;

  std::vector<double> forces_eva(n_atoms * 3u, 0.0);
  lammps_gather_atoms(session->lmp, "f", 1, 3, forces_eva.data());
  for (size_t i = 0; i < forces_eva.size(); ++i)
    forces_h_bohr[i] = forces_eva[i] * kAngstromPerBohr / kEvPerHartree;

  const double pe_ev = lammps_get_thermo(session->lmp, "pe");
  *energy_h = pe_ev / kEvPerHartree;

  const double sxx = lammps_get_thermo(session->lmp, "pxx");
  const double syy = lammps_get_thermo(session->lmp, "pyy");
  const double szz = lammps_get_thermo(session->lmp, "pzz");
  const double sxy = lammps_get_thermo(session->lmp, "pxy");
  const double sxz = lammps_get_thermo(session->lmp, "pxz");
  const double syz = lammps_get_thermo(session->lmp, "pyz");
  const double evA3_to_hbohr3 =
      (1.0 / kEvPerHartree) * std::pow(kAngstromPerBohr, 3.0);
  const double c = -(1.0 / kBarPerEvA3) * evA3_to_hbohr3;
  stress_h_bohr3[0] = sxx * c;
  stress_h_bohr3[1] = sxy * c;
  stress_h_bohr3[2] = sxz * c;
  stress_h_bohr3[3] = sxy * c;
  stress_h_bohr3[4] = syy * c;
  stress_h_bohr3[5] = syz * c;
  stress_h_bohr3[6] = sxz * c;
  stress_h_bohr3[7] = syz * c;
  stress_h_bohr3[8] = szz * c;
  return 0;
}

int reject_common_overlay(CommonMethodSpec_ptr common) {
  if (common.p.type == CAPN_NULL)
    return 0;
  struct CommonMethodSpec view;
  read_CommonMethodSpec(&view, common);
  capn_resolve(&view.xcFunctionals);
  capn_resolve(&view.kMesh.p);
  int has_any = view.xcFunctionals.len > 0 || view.basisSet.len > 0 ||
                view.planewaveCutoffEv > 0.0 || view.charge != 0 ||
                view.spinMultiplicity > 0 || view.scfEnergyToleranceEv > 0.0 ||
                view.scfMaxIterations > 0 || view.kMesh.p.len > 0 ||
                view.vanDerWaalsMethod.len > 0 ||
                view.relativityMethod.len > 0 || view.vanDerWaalsS6 > 0.0;
  capn_resolve(&view.smearing.p);
  if (view.smearing.p.type == CAPN_STRUCT) {
    struct CommonMethodSpec_Smearing smear;
    read_CommonMethodSpec_Smearing(&smear, view.smearing);
    if (smear.kind != CommonMethodSpec_Smearing_Kind_none)
      has_any = 1;
  }
  if (has_any) {
    store_error("common overlay is unsupported by lammpc");
    return -1;
  }
  return 0;
}

LammpcSession *session_from_params_ptr(LammpsParams_ptr params) {
  auto *session = new LammpcSession();
  if (decode_lammps_params(params, &session->params) != 0) {
    delete session;
    return nullptr;
  }
  clear_error();
  return session;
}

LammpcSession *session_from_config_bytes(const void *config_capnp,
                                         size_t config_capnp_size_bytes) {
  struct capn arena;
  capn_ptr root;
  if (!config_capnp || config_capnp_size_bytes == 0)
    return nullptr;
  memset(&arena, 0, sizeof(arena));
  if (capn_init_mem(&arena, (const uint8_t *)config_capnp,
                    config_capnp_size_bytes, 0) != 0) {
    store_error("invalid PotentialConfig message");
    return nullptr;
  }
  root = capn_getp(capn_root(&arena), 0, 1);
  if (root.type != CAPN_STRUCT) {
    capn_free(&arena);
    store_error("invalid PotentialConfig message");
    return nullptr;
  }
  PotentialConfig_ptr config;
  config.p = root;
  struct PotentialConfig view;
  read_PotentialConfig(&view, config);
  if (reject_common_overlay(view.common) != 0) {
    capn_free(&arena);
    return nullptr;
  }
  if (view.which != PotentialConfig_lammps) {
    capn_free(&arena);
    store_error("PotentialConfig does not carry the lammps arm");
    return nullptr;
  }
  LammpcSession *session = session_from_params_ptr(view.lammps);
  capn_free(&arena);
  return session;
}

extern "C" {

const char *lammpc_last_error(void) { return g_last_error; }

int lammpc_set_params(const void *, size_t) {
  store_error("lammpc: use sessions; global params are not supported");
  return -1;
}

int lammpc_configure(const void *config_capnp, size_t config_capnp_size_bytes) {
  LammpcSession *session =
      session_from_config_bytes(config_capnp, config_capnp_size_bytes);
  if (!session)
    return -1;
  lammpc_session_destroy(session);
  return 0;
}

LammpcSession *lammpc_session_create(const void *params_capnp,
                                     size_t params_capnp_size_bytes) {
  struct capn arena;
  LammpsParams_ptr params;
  if (lammpc_params_root(params_capnp, params_capnp_size_bytes, &arena,
                         &params) != 0) {
    store_error("invalid LammpsParams message");
    return nullptr;
  }
  LammpcSession *session = session_from_params_ptr(params);
  lammpc_params_release(&arena);
  return session;
}

int lammpc_session_set_params(LammpcSession *session, const void *params_capnp,
                              size_t params_capnp_size_bytes) {
  if (!session)
    return -1;
  struct capn arena;
  LammpsParams_ptr params;
  if (lammpc_params_root(params_capnp, params_capnp_size_bytes, &arena,
                         &params) != 0) {
    store_error("invalid LammpsParams message");
    return -1;
  }
  DecodedParams decoded;
  int rc = decode_lammps_params(params, &decoded);
  lammpc_params_release(&arena);
  if (rc != 0)
    return -1;
  session->params = decoded;
  destroy_instance(session);
  clear_error();
  return 0;
}

LammpcSession *lammpc_session_create_from_config(const void *config_capnp,
                                                 size_t config_capnp_size_bytes) {
  return session_from_config_bytes(config_capnp, config_capnp_size_bytes);
}

int lammpc_session_configure(LammpcSession *session, const void *config_capnp,
                             size_t config_capnp_size_bytes) {
  if (!session)
    return -1;
  LammpcSession *fresh =
      session_from_config_bytes(config_capnp, config_capnp_size_bytes);
  if (!fresh)
    return -1;
  session->params = fresh->params;
  destroy_instance(session);
  lammpc_session_destroy(fresh);
  clear_error();
  return 0;
}

void lammpc_session_destroy(LammpcSession *session) {
  if (!session)
    return;
  destroy_instance(session);
  delete session;
}

LammpcResult lammpc_session_calculate_result(
    LammpcSession *session, const void *force_input_capnp,
    size_t force_input_capnp_size_bytes, void *potential_result_capnp,
    size_t potential_result_capnp_capacity_bytes,
    size_t *potential_result_capnp_size_bytes) {
  if (!session || !potential_result_capnp_size_bytes)
    return fail_result("invalid arguments");
  *potential_result_capnp_size_bytes = 0;

  struct capn arena;
  ForceInput_ptr force_input;
  if (lammpc_force_input_root(force_input_capnp, force_input_capnp_size_bytes,
                              &arena, &force_input) != 0)
    return fail_result("invalid ForceInput message");

  size_t n_atoms = 0;
  int has_cell = 0;
  double energy_factor = 1.0;
  double force_factor = 1.0;
  if (lammpc_force_input_atom_count(force_input, &n_atoms, &has_cell) != 0 ||
      lammpc_force_input_result_factors(force_input, &energy_factor,
                                        &force_factor) != 0) {
    lammpc_params_release(&arena);
    return fail_result("invalid ForceInput geometry or units");
  }
  size_t required = lammpc_potential_result_flat_size(n_atoms * 3u);
  *potential_result_capnp_size_bytes = required;
  if (!potential_result_capnp ||
      potential_result_capnp_capacity_bytes < required) {
    lammpc_params_release(&arena);
    return fail_result("PotentialResult buffer too small");
  }

  std::vector<double> positions_ang(n_atoms * 3u);
  std::vector<int> atomic_numbers(n_atoms);
  double cell_ang[9];
  if (lammpc_force_input_copy_geometry(force_input, positions_ang.data(),
                                       atomic_numbers.data(), n_atoms,
                                       cell_ang, &has_cell) != 0) {
    lammpc_params_release(&arena);
    return fail_result("invalid ForceInput geometry");
  }
  lammpc_params_release(&arena);

  if (!geometry_matches(session, n_atoms, cell_ang, has_cell)) {
    if (build_instance(session, n_atoms, positions_ang.data(),
                       atomic_numbers.data(), cell_ang, has_cell) != 0)
      return fail_result(g_last_error);
  }

  double energy_h = 0.0;
  std::vector<double> forces_h_bohr(n_atoms * 3u);
  double stress_h_bohr3[9] = {0};
  if (evaluate(session, n_atoms, positions_ang.data(), &energy_h,
               forces_h_bohr.data(), stress_h_bohr3) != 0)
    return fail_result(g_last_error);

  std::vector<double> forces_out(n_atoms * 3u);
  for (size_t i = 0; i < forces_out.size(); ++i)
    forces_out[i] = forces_h_bohr[i] * force_factor;
  const double length_factor =
      energy_factor != 0.0
          ? force_factor * kAngstromPerBohr / energy_factor
          : 1.0;
  const double stress_factor =
      energy_factor * std::pow(length_factor / kAngstromPerBohr, 3.0);
  double stress_out[9];
  for (int i = 0; i < 9; ++i)
    stress_out[i] = stress_h_bohr3[i] * stress_factor;

  if (lammpc_potential_result_write(energy_h * energy_factor, forces_out.data(),
                                    forces_out.size(), stress_out,
                                    potential_result_capnp,
                                    potential_result_capnp_capacity_bytes,
                                    potential_result_capnp_size_bytes) != 0)
    return fail_result("PotentialResult write failed");

  LammpcResult r;
  r.ok = 1;
  r.energy_h = energy_h;
  std::snprintf(r.message, sizeof(r.message), "ok");
  return r;
}

LammpcResult lammpc_calculate_result(const void *params_capnp,
                                     size_t params_capnp_size_bytes,
                                     const void *force_input_capnp,
                                     size_t force_input_capnp_size_bytes,
                                     void *potential_result_capnp,
                                     size_t potential_result_capnp_capacity_bytes,
                                     size_t *potential_result_capnp_size_bytes) {
  LammpcSession *session =
      lammpc_session_create(params_capnp, params_capnp_size_bytes);
  if (!session)
    return fail_result(g_last_error);
  LammpcResult r = lammpc_session_calculate_result(
      session, force_input_capnp, force_input_capnp_size_bytes,
      potential_result_capnp, potential_result_capnp_capacity_bytes,
      potential_result_capnp_size_bytes);
  lammpc_session_destroy(session);
  return r;
}

LammpcResult lammpc_calculate_result_from_config(
    const void *config_capnp, size_t config_capnp_size_bytes,
    const void *force_input_capnp, size_t force_input_capnp_size_bytes,
    void *potential_result_capnp,
    size_t potential_result_capnp_capacity_bytes,
    size_t *potential_result_capnp_size_bytes) {
  LammpcSession *session =
      lammpc_session_create_from_config(config_capnp, config_capnp_size_bytes);
  if (!session)
    return fail_result(g_last_error);
  LammpcResult r = lammpc_session_calculate_result(
      session, force_input_capnp, force_input_capnp_size_bytes,
      potential_result_capnp, potential_result_capnp_capacity_bytes,
      potential_result_capnp_size_bytes);
  lammpc_session_destroy(session);
  return r;
}

size_t lammpc_potential_result_size_for_force_input(
    const void *force_input_capnp, size_t force_input_capnp_size_bytes) {
  struct capn arena;
  ForceInput_ptr force_input;
  if (lammpc_force_input_root(force_input_capnp, force_input_capnp_size_bytes,
                              &arena, &force_input) != 0)
    return 0;
  size_t n_atoms = 0;
  int has_cell = 0;
  int rc = lammpc_force_input_atom_count(force_input, &n_atoms, &has_cell);
  lammpc_params_release(&arena);
  if (rc != 0)
    return 0;
  return lammpc_potential_result_flat_size(n_atoms * 3u);
}

LammpcResult lammpc_energy_gradient(int n_atoms, const double *positions_ang,
                                    const int *atomic_numbers,
                                    const void *cell_ang, size_t cell_len,
                                    double *gradient_h_bohr) {
  LammpcResult r =
      lammpc_energy_forces(n_atoms, positions_ang, atomic_numbers, cell_ang,
                           cell_len, gradient_h_bohr);
  if (r.ok && gradient_h_bohr)
    for (int i = 0; i < 3 * n_atoms; ++i)
      gradient_h_bohr[i] = -gradient_h_bohr[i];
  return r;
}

LammpcResult lammpc_energy(int n_atoms, const double *positions_ang,
                           const int *atomic_numbers, const void *cell_ang,
                           size_t cell_len) {
  return lammpc_energy_forces(n_atoms, positions_ang, atomic_numbers, cell_ang,
                              cell_len, nullptr);
}

LammpcResult lammpc_energy_forces(int, const double *, const int *, const void *,
                                  size_t, double *) {
  return fail_result("lammpc: raw-array calls need a session");
}

const char *lammpc_version(void) { return "lammpc/api"; }
int lammpc_abi_version(void) { return LAMMPC_ABI_VERSION; }
int lammpc_available(void) { return 1; }
void lammpc_finalize(void) {}

} // extern "C"
