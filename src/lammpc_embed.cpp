/**
 * lammpc embed: direct in-process calls to the actual liblammps functions.
 *
 * Setup replicates what LAMMPS' input commands do by calling the same
 * member functions they call (units via Update::set_units, pair creation via
 * Force::create_pair, coefficients via Pair::settings/coeff, box fields plus
 * Domain::set_initial_box/set_global_box/set_local_box, atoms via
 * AtomVec::create_atom). The step loop writes positions straight into
 * Atom::x, calls Pair::compute, and reads Atom::f and the pair energy
 * accumulators. No command-string parsing, no input scripts, no file
 * traffic anywhere.
 */
#include "lammpc.h"
#include "lammpc_params.h"

#include <mpi.h>

#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "domain.h"
#include "force.h"
#include "integrate.h"
#include "lammps.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "pair.h"
#include "update.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#ifndef LAMMPC_VERSION_STRING
#define LAMMPC_VERSION_STRING "unknown"
#endif

namespace {

thread_local char g_last_error[512] = "";

void store_error(const char *message) {
  std::snprintf(g_last_error, sizeof(g_last_error), "%s", message);
}

void clear_error() { g_last_error[0] = '\0'; }

LammpcResult fail_result(const char *message) {
  LammpcResult r;
  r.ok = 0;
  r.energy_h = 0.0;
  std::snprintf(r.message, sizeof(r.message), "%s", message);
  return r;
}

/* Tokenize one of OUR schema text fields (pairStyle, pairCoeffs entries,
 * boundary) into argv form for the LAMMPS member functions. This tokenizes
 * lammpc's own wire fields; LAMMPS never parses command lines here. */
std::vector<std::string> split_tokens(const char *text, int len) {
  std::vector<std::string> out;
  std::string current;
  for (int i = 0; i < len; ++i) {
    char c = text[i];
    if (c == ' ' || c == '\t') {
      if (!current.empty()) {
        out.push_back(current);
        current.clear();
      }
    } else {
      current.push_back(c);
    }
  }
  if (!current.empty())
    out.push_back(current);
  return out;
}

std::vector<char *> argv_view(std::vector<std::string> &tokens) {
  std::vector<char *> argv;
  argv.reserve(tokens.size());
  for (auto &token : tokens)
    argv.push_back(token.data());
  return argv;
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
  bool newton_pair = true;
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
  out->newton_pair = view.newtonPair != 0;

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

/* Standard atomic weights (amu), Z-indexed, H..Og subset used when
 * LammpsParams.masses is empty. Zero marks elements callers must supply. */
const double kStandardMassAmu[119] = {
    0,      1.008,  4.0026, 6.94,   9.0122, 10.81,  12.011, 14.007, 15.999,
    18.998, 20.180, 22.990, 24.305, 26.982, 28.085, 30.974, 32.06,  35.45,
    39.948, 39.098, 40.078, 44.956, 47.867, 50.942, 51.996, 54.938, 55.845,
    58.933, 58.693, 63.546, 65.38,  69.723, 72.630, 74.922, 78.971, 79.904,
    83.798, 85.468, 87.62,  88.906, 91.224, 92.906, 95.95,  97.0,   101.07,
    102.91, 106.42, 107.87, 112.41, 114.82, 118.71, 121.76, 127.60, 126.90,
    131.29, 132.91, 137.33, 138.91, 140.12, 140.91, 144.24, 145.0,  150.36,
    151.96, 157.25, 158.93, 162.50, 164.93, 167.26, 168.93, 173.05, 174.97,
    178.49, 180.95, 183.84, 186.21, 190.23, 192.22, 195.08, 196.97, 200.59,
    204.38, 207.2,  208.98, 209.0,  210.0,  222.0,  223.0,  226.0,  227.0,
    232.04, 231.04, 238.03, 237.0,  244.0,  243.0,  247.0,  247.0,  251.0,
    252.0,  257.0,  258.0,  259.0,  262.0,  267.0,  270.0,  269.0,  270.0,
    270.0,  278.0,  281.0,  281.0,  285.0,  286.0,  289.0,  289.0,  293.0,
    293.0,  294.0};

constexpr double kEvPerHartree = 27.211386245988;
constexpr double kAngstromPerBohr = 0.529177210903;

} // namespace

struct LammpcSession {
  LAMMPS_NS::LAMMPS *lmp = nullptr;
  DecodedParams params;
  /* Geometry signature the live instance was built for. */
  size_t n_atoms = 0;
  double cell_ang[9] = {0};
  bool box_ready = false;
};

namespace {

void destroy_instance(LammpcSession *session) {
  if (session->lmp) {
    delete session->lmp;
    session->lmp = nullptr;
  }
  session->box_ready = false;
}

int ensure_mpi() {
  int initialized = 0;
  MPI_Initialized(&initialized);
  if (!initialized) {
    int argc = 0;
    char **argv = nullptr;
    if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
      store_error("MPI_Init failed");
      return -1;
    }
  }
  return 0;
}

int type_for_atomic_number(const DecodedParams &params, int z) {
  for (size_t t = 0; t < params.type_to_z.size(); ++t) {
    if (params.type_to_z[t] == z)
      return (int)t + 1;
  }
  return -1;
}

/* Build one LAMMPS instance for this geometry through the same member
 * functions the create_box / create_atoms / pair_style / pair_coeff
 * commands call. */
int build_instance(LammpcSession *session, size_t n_atoms,
                   const double *positions_ang, const int *atomic_numbers,
                   const double *cell_ang, int has_cell) {
  using namespace LAMMPS_NS;
  destroy_instance(session);
  if (ensure_mpi() != 0)
    return -1;

  const DecodedParams &p = session->params;
  const char *ctor_args[] = {"lammpc", "-log", "none", "-screen", "none",
                             "-nocite"};
  LAMMPS *lmp = nullptr;
  try {
    lmp = new LAMMPS(6, const_cast<char **>(ctor_args), MPI_COMM_WORLD);

    lmp->update->set_units(p.units_style.c_str());

    {
      auto boundary_tokens = split_tokens(p.boundary.c_str(),
                                          (int)p.boundary.size());
      if (boundary_tokens.size() != 3) {
        store_error("LammpsParams.boundary needs 3 tokens");
        delete lmp;
        return -1;
      }
      auto argv = argv_view(boundary_tokens);
      lmp->domain->set_boundary(3, argv.data(), 0);
    }

    lmp->atom->create_avec(p.atom_style.c_str(), 0, nullptr, 1);
    lmp->atom->ntypes = (int)p.type_to_z.size();

    /* Box: fields plus the same call sequence CreateBox::command runs. */
    double lo[3] = {0, 0, 0};
    double hi[3] = {0, 0, 0};
    if (has_cell) {
      for (int d = 0; d < 3; ++d) {
        lo[d] = 0.0;
        hi[d] = cell_ang[4 * d]; /* diagonal of the row-major 3x3 cell */
        if (!(hi[d] > 0.0)) {
          store_error("only diagonal positive cells are supported");
          delete lmp;
          return -1;
        }
      }
      if (cell_ang[1] != 0.0 || cell_ang[2] != 0.0 || cell_ang[3] != 0.0 ||
          cell_ang[5] != 0.0 || cell_ang[6] != 0.0 || cell_ang[7] != 0.0) {
        store_error("only diagonal cells are supported");
        delete lmp;
        return -1;
      }
    } else {
      /* Shrink-wrapped styles come from boundary; still need an initial
       * extent covering the atoms. */
      double pad = 50.0;
      double mins[3] = {1e300, 1e300, 1e300};
      double maxs[3] = {-1e300, -1e300, -1e300};
      for (size_t i = 0; i < n_atoms; ++i) {
        for (int d = 0; d < 3; ++d) {
          double v = positions_ang[3 * i + d];
          if (v < mins[d])
            mins[d] = v;
          if (v > maxs[d])
            maxs[d] = v;
        }
      }
      for (int d = 0; d < 3; ++d) {
        lo[d] = mins[d] - pad;
        hi[d] = maxs[d] + pad;
      }
    }
    for (int d = 0; d < 3; ++d) {
      lmp->domain->boxlo[d] = lo[d];
      lmp->domain->boxhi[d] = hi[d];
    }
    lmp->domain->box_exist = 1;
    lmp->domain->set_initial_box();
    lmp->domain->set_global_box();
    lmp->comm->set_proc_grid();
    lmp->domain->set_local_box();

    /* Masses per type. */
    for (int t = 1; t <= lmp->atom->ntypes; ++t) {
      double mass = 0.0;
      if ((size_t)t <= p.masses.size())
        mass = p.masses[(size_t)t - 1];
      if (mass <= 0.0) {
        int z = p.type_to_z[(size_t)t - 1];
        if (z > 0 && z < 119)
          mass = kStandardMassAmu[z];
      }
      if (mass <= 0.0) {
        store_error("no mass for a LAMMPS type (set LammpsParams.masses)");
        delete lmp;
        return -1;
      }
      lmp->atom->set_mass("lammpc", 0, t, mass);
    }

    /* Atoms in caller order (tags 1..N follow creation order). */
    for (size_t i = 0; i < n_atoms; ++i) {
      int type = type_for_atomic_number(p, atomic_numbers[i]);
      if (type < 1) {
        store_error("ForceInput atomic number missing from "
                    "LammpsParams.typeToAtomicNumber");
        delete lmp;
        return -1;
      }
      double coord[3] = {positions_ang[3 * i + 0], positions_ang[3 * i + 1],
                         positions_ang[3 * i + 2]};
      lmp->atom->avec->create_atom(type, coord);
    }
    lmp->atom->natoms = (double)n_atoms;
    lmp->atom->tag_extend();

    lmp->force->newton = lmp->force->newton_pair = lmp->force->newton_bond =
        p.newton_pair ? 1 : 0;

    /* Pair style + coefficients through the actual member functions. */
    {
      auto style_tokens =
          split_tokens(p.pair_style.c_str(), (int)p.pair_style.size());
      if (style_tokens.empty()) {
        store_error("empty pair style");
        delete lmp;
        return -1;
      }
      std::string style_name = style_tokens.front();
      lmp->force->create_pair(style_name.c_str(), 1);
      std::vector<std::string> setting_tokens(style_tokens.begin() + 1,
                                              style_tokens.end());
      auto argv = argv_view(setting_tokens);
      lmp->force->pair->settings((int)setting_tokens.size(), argv.data());
    }
    for (const auto &coeff_line : session->params.pair_coeffs) {
      auto coeff_tokens =
          split_tokens(coeff_line.c_str(), (int)coeff_line.size());
      auto argv = argv_view(coeff_tokens);
      lmp->force->pair->coeff((int)coeff_tokens.size(), argv.data());
    }

    /* Init + neighbor build + first force evaluation. */
    lmp->update->whichflag = 1;
    lmp->init();
    lmp->update->integrate->setup(0);
    lmp->update->whichflag = 0;
  } catch (std::exception &ex) {
    store_error(ex.what());
    delete lmp;
    return -1;
  }

  session->lmp = lmp;
  session->n_atoms = n_atoms;
  if (has_cell)
    std::memcpy(session->cell_ang, cell_ang, sizeof(session->cell_ang));
  else
    std::memset(session->cell_ang, 0, sizeof(session->cell_ang));
  session->box_ready = true;
  return 0;
}

bool geometry_matches(const LammpcSession *session, size_t n_atoms,
                      const double *cell_ang, int has_cell) {
  if (!session->box_ready || session->n_atoms != n_atoms)
    return false;
  for (int i = 0; i < 9; ++i) {
    double want = has_cell ? cell_ang[i] : 0.0;
    if (session->cell_ang[i] != want)
      return false;
  }
  return true;
}

/* One energy+forces evaluation on a live instance: positions in, direct
 * pair compute, energy + per-tag forces out (Hartree / Hartree/Bohr). */
int evaluate(LammpcSession *session, size_t n_atoms,
             const double *positions_ang, double *energy_h,
             double *forces_h_bohr) {
  using namespace LAMMPS_NS;
  LAMMPS *lmp = session->lmp;
  if (!lmp || lmp->atom->nlocal != (int)n_atoms) {
    store_error("session instance does not match the step geometry");
    return -1;
  }

  Atom *atom = lmp->atom;
  double **x = atom->x;
  /* Positions per original tag: creation order gave tag i+1 to caller
   * index i, and tags survive any internal reordering. */
  for (int i = 0; i < atom->nlocal; ++i) {
    int caller = (int)atom->tag[i] - 1;
    if (caller < 0 || caller >= (int)n_atoms) {
      store_error("unexpected atom tag in the live instance");
      return -1;
    }
    x[i][0] = positions_ang[3 * caller + 0];
    x[i][1] = positions_ang[3 * caller + 1];
    x[i][2] = positions_ang[3 * caller + 2];
  }

  try {
    /* Re-neighbor through the same path a timestep uses. */
    lmp->update->whichflag = 1;
    lmp->update->integrate->setup_minimal(1);
    lmp->update->whichflag = 0;

    /* Direct force pass with the global energy accumulators enabled. */
    int nall = atom->nlocal + atom->nghost;
    if (nall > 0)
      std::memset(&atom->f[0][0], 0, sizeof(double) * 3 * (size_t)nall);
    lmp->force->pair->compute(1, 0);
    if (lmp->force->newton_pair)
      lmp->comm->reverse_comm();
  } catch (std::exception &ex) {
    store_error(ex.what());
    return -1;
  }

  double energy_ev =
      lmp->force->pair->eng_vdwl + lmp->force->pair->eng_coul;
  *energy_h = energy_ev / kEvPerHartree;
  if (forces_h_bohr) {
    double **f = lmp->atom->f;
    for (int i = 0; i < lmp->atom->nlocal; ++i) {
      int caller = (int)lmp->atom->tag[i] - 1;
      for (int d = 0; d < 3; ++d)
        forces_h_bohr[3 * caller + d] =
            f[i][d] * kAngstromPerBohr / kEvPerHartree;
    }
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

/* Reject any set CommonMethodSpec field: a classical force field lowers
 * none of them, and silence would misconfigure the caller. */
int reject_common_overlay(CommonMethodSpec_ptr common) {
  if (common.p.type == CAPN_NULL)
    return 0;
  struct CommonMethodSpec view;
  read_CommonMethodSpec(&view, common);
  capn_resolve(&view.xcFunctionals);
  capn_resolve(&view.kMesh.p);
  int has_any = view.xcFunctionals.len > 0 || view.basisSet.len > 0 ||
                view.planewaveCutoffEv > 0.0 || view.charge != 0 ||
                view.spinMultiplicity > 0 ||
                view.scfEnergyToleranceEv > 0.0 ||
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
    store_error("common overlay: no CommonMethodSpec field has a LAMMPS "
                "lowering; configure via the lammps arm");
    return -1;
  }
  return 0;
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

} // namespace

extern "C" {

const char *lammpc_last_error(void) { return g_last_error; }

int lammpc_set_params(const void *params_capnp,
                      size_t params_capnp_size_bytes) {
  (void)params_capnp;
  (void)params_capnp_size_bytes;
  store_error("lammpc: use sessions; global params are not supported");
  return -1;
}

int lammpc_configure(const void *config_capnp,
                     size_t config_capnp_size_bytes) {
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

LammpcSession *lammpc_session_create_from_config(
    const void *config_capnp, size_t config_capnp_size_bytes) {
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
  if (lammpc_force_input_root(force_input_capnp,
                              force_input_capnp_size_bytes, &arena,
                              &force_input) != 0)
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
  if (evaluate(session, n_atoms, positions_ang.data(), &energy_h,
               forces_h_bohr.data()) != 0)
    return fail_result(g_last_error);

  std::vector<double> forces_out(n_atoms * 3u);
  for (size_t i = 0; i < forces_out.size(); ++i)
    forces_out[i] = forces_h_bohr[i] * force_factor;
  if (lammpc_potential_result_write(energy_h * energy_factor,
                                    forces_out.data(), forces_out.size(),
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
  if (lammpc_force_input_root(force_input_capnp,
                              force_input_capnp_size_bytes, &arena,
                              &force_input) != 0)
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
  if (r.ok && gradient_h_bohr) {
    for (int i = 0; i < 3 * n_atoms; ++i)
      gradient_h_bohr[i] = -gradient_h_bohr[i];
  }
  return r;
}

LammpcResult lammpc_energy(int n_atoms, const double *positions_ang,
                           const int *atomic_numbers, const void *cell_ang,
                           size_t cell_len) {
  return lammpc_energy_forces(n_atoms, positions_ang, atomic_numbers,
                              cell_ang, cell_len, nullptr);
}

LammpcResult lammpc_energy_forces(int n_atoms, const double *positions_ang,
                                  const int *atomic_numbers,
                                  const void *cell_ang, size_t cell_len,
                                  double *forces_h_bohr) {
  (void)n_atoms;
  (void)positions_ang;
  (void)atomic_numbers;
  (void)cell_ang;
  (void)cell_len;
  (void)forces_h_bohr;
  return fail_result(
      "lammpc: raw-array calls need a session (no global params)");
}

const char *lammpc_version(void) { return "lammpc/" LAMMPC_VERSION_STRING; }

int lammpc_abi_version(void) { return LAMMPC_ABI_VERSION; }

int lammpc_available(void) { return 1; }

void lammpc_finalize(void) {}

} // extern "C"
