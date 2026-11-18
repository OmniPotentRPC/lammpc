#include "lammpc.h"
#include "lammpc_params.h"
#include "lammpc_boundary_embed.h"

#include "atom.h"
#include "domain.h"
#include "lammps.h"
#include "library.h"
#include "minimage.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <sstream>
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
/* LAMMPS metal-unit conversion: mvv2e * amu * (angstrom/ps)^2 = eV. */
constexpr double kMassVelocityToEv = 1.0364269656262175e-4;

struct HardBoundaryHook {
  LammpcHardBoundaryFn fn = nullptr;
  void *user = nullptr;
};

struct LammpcSession {
  void *lmp = nullptr;
  DecodedParams params;
  size_t n_atoms = 0;
  int has_cell = 0;
  double cell_ang[9] = {0};
  std::vector<int> atom_types;
  bool md_ready = false;
  LammpcMdConfig md_config = {};
  std::vector<double> md_point;
  std::vector<double> md_normal;
  std::vector<double> md_masses;
  std::vector<double> md_positions;
  std::vector<double> md_velocities;
  std::vector<double> md_forces;
  std::vector<double> md_removed_velocities;
  double md_constraint_force = 0.0;
  size_t md_force_evaluations = 0;
  size_t md_steps = 0;
  uint64_t md_rng[4] = {0, 0, 0, 0};
  bool nve_fix_installed = false;
  std::vector<double> nve_positions;
  std::vector<double> nve_velocities;
  LammpcNveObservation nve_observation = {};
  HardBoundaryHook hard_boundary;
};

struct LammpsBounds {
  double xlo = 0.0;
  double xhi = 0.0;
  double ylo = 0.0;
  double yhi = 0.0;
  double zlo = 0.0;
  double zhi = 0.0;
  double xy = 0.0;
  double xz = 0.0;
  double yz = 0.0;
  bool triclinic = false;
};

std::string format_double17(double value) {
  char buffer[64] = {0};
  std::snprintf(buffer, sizeof(buffer), "%.17g", value);
  return std::string(buffer);
}

std::string minimage_last_error_or(const char *fallback) {
  const char *msg = mi_last_error();
  if (msg && msg[0] != '\0')
    return std::string(msg);
  return std::string(fallback);
}

bool boundary_is_fully_periodic(const std::string &boundary) {
  std::istringstream tokens(boundary);
  std::string bx, by, bz;
  if (!(tokens >> bx >> by >> bz))
    return false;
  auto is_periodic = [](const std::string &token) {
    return !token.empty() && (token[0] == 'p' || token[0] == 'P');
  };
  return is_periodic(bx) && is_periodic(by) && is_periodic(bz);
}

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
  session->md_ready = false;
  session->nve_fix_installed = false;
  session->nve_positions.clear();
  session->nve_velocities.clear();
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
  LammpsBounds bounds;
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
    bounds.xlo = 0.0;
    bounds.xhi = cell_ang[0];
    bounds.ylo = 0.0;
    bounds.yhi = cell_ang[4];
    bounds.zlo = 0.0;
    bounds.zhi = cell_ang[8];
    bounds.xy = cell_ang[3];
    bounds.xz = cell_ang[6];
    bounds.yz = cell_ang[7];
    bounds.triclinic = true;
    region_cmd = "region box prism " + format_double17(bounds.xlo) + " " +
                 format_double17(bounds.xhi) + " " +
                 format_double17(bounds.ylo) + " " +
                 format_double17(bounds.yhi) + " " +
                 format_double17(bounds.zlo) + " " +
                 format_double17(bounds.zhi) + " " +
                 format_double17(bounds.xy) + " " +
                 format_double17(bounds.xz) + " " +
                 format_double17(bounds.yz) + " units box";
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
    bounds.xlo = mins[0] - pad;
    bounds.xhi = maxs[0] + pad;
    bounds.ylo = mins[1] - pad;
    bounds.yhi = maxs[1] + pad;
    bounds.zlo = mins[2] - pad;
    bounds.zhi = maxs[2] + pad;
    bounds.xy = 0.0;
    bounds.xz = 0.0;
    bounds.yz = 0.0;
    bounds.triclinic = false;
    region_cmd = "region box block " + format_double17(bounds.xlo) + " " +
                 format_double17(bounds.xhi) + " " +
                 format_double17(bounds.ylo) + " " +
                 format_double17(bounds.yhi) + " " +
                 format_double17(bounds.zlo) + " " +
                 format_double17(bounds.zhi) + " units box";
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
                         format_double17(mass)) != 0) {
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
  }

  if (n_atoms > (size_t)std::numeric_limits<int>::max()) {
    store_error("atom count exceeds lammps_create_atoms int range");
    lammps_close(lmp);
    return -1;
  }

  std::vector<double> creation_positions(positions_ang,
                                         positions_ang + n_atoms * 3u);
  if (n_atoms > 0 && boundary_is_fully_periodic(p.boundary)) {
    mi_cell minimage_cell{};
    if (mi_cell_from_lammps(bounds.xlo, bounds.xhi, bounds.ylo, bounds.yhi,
                            bounds.zlo, bounds.zhi, bounds.xy, bounds.xz,
                            bounds.yz, &minimage_cell) != 0) {
      store_error("minimage cell construction failed: " +
                  minimage_last_error_or("unknown minimage error"));
      lammps_close(lmp);
      return -1;
    }
    const double center[3] = {
        minimage_cell.ox + 0.5 * (minimage_cell.ax + minimage_cell.bx +
                                  minimage_cell.cx),
        minimage_cell.oy + 0.5 * (minimage_cell.ay + minimage_cell.by +
                                  minimage_cell.cy),
        minimage_cell.oz + 0.5 * (minimage_cell.az + minimage_cell.bz +
                                  minimage_cell.cz),
    };
    std::vector<double> centered(n_atoms * 3u, 0.0);
    std::vector<double> wrapped(n_atoms * 3u, 0.0);
    for (size_t i = 0; i < n_atoms; ++i) {
      centered[3 * i + 0] = creation_positions[3 * i + 0] - center[0];
      centered[3 * i + 1] = creation_positions[3 * i + 1] - center[1];
      centered[3 * i + 2] = creation_positions[3 * i + 2] - center[2];
    }
    if (mi_wrap_many(&minimage_cell, centered.data(), n_atoms, wrapped.data()) !=
        0) {
      store_error("minimage coordinate wrapping failed: " +
                  minimage_last_error_or("unknown minimage error"));
      lammps_close(lmp);
      return -1;
    }
    for (size_t i = 0; i < n_atoms; ++i) {
      for (int d = 0; d < 3; ++d) {
        const size_t k = 3 * i + (size_t)d;
        const double candidate = center[d] + wrapped[k];
        const double scale = std::max(
            1.0, std::max(std::fabs(creation_positions[k]), std::fabs(candidate)));
        const double tolerance =
            128.0 * std::numeric_limits<double>::epsilon() * scale;
        if (std::fabs(creation_positions[k] - candidate) > tolerance)
          creation_positions[k] = candidate;
      }
    }
  }

#if defined(LAMMPS_BIGBIG)
  std::vector<int64_t> image_flags(
      n_atoms, lammps_encode_image_flags(0, 0, 0));
#else
  std::vector<int> image_flags(n_atoms, lammps_encode_image_flags(0, 0, 0));
#endif
  int created = 0;
  if (n_atoms > 0) {
    created = lammps_create_atoms(
        lmp, (int)n_atoms, nullptr, session->atom_types.data(),
        creation_positions.data(), nullptr, image_flags.data(), 1);
  }
  if (created < 0 || lammps_has_error(lmp)) {
    char buffer[512] = {0};
    if (lammps_get_last_error_message(lmp, buffer, (int)sizeof(buffer)) > 0)
      store_error(buffer);
    else
      store_error("lammps_create_atoms failed");
    lammps_close(lmp);
    return -1;
  }
  const long long present = (long long)std::llround(lammps_get_natoms(lmp));
  if (created != (int)n_atoms || present != (long long)n_atoms) {
    char mismatch[256] = {0};
    std::snprintf(mismatch, sizeof(mismatch),
                  "atom count mismatch after creation: requested=%zu "
                  "created=%d present=%lld",
                  n_atoms, created, present);
    store_error(mismatch);
    lammps_close(lmp);
    return -1;
  }
  if (n_atoms > 0) {
    lammps_scatter_atoms(lmp, "x", 1, 3, creation_positions.data());
    if (lammps_has_error(lmp)) {
      char buffer[512] = {0};
      if (lammps_get_last_error_message(lmp, buffer, (int)sizeof(buffer)) > 0)
        store_error(buffer);
      else
        store_error("lammps_scatter_atoms after creation failed");
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
  if (run_cmd(lmp, "thermo_style custom step pe ke pxx pyy pzz pxy pxz pyz") != 0 ||
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

uint64_t splitmix64(uint64_t *state) {
  *state += UINT64_C(0x9e3779b97f4a7c15);
  uint64_t z = *state;
  z = (z ^ (z >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
  z = (z ^ (z >> 27)) * UINT64_C(0x94d049bb133111eb);
  return z ^ (z >> 31);
}

void md_seed(LammpcSession *session, unsigned long long value) {
  uint64_t seed = (uint64_t)value;
  for (uint64_t &slot : session->md_rng)
    slot = splitmix64(&seed);
}

uint64_t md_random_u64(LammpcSession *session) {
  uint64_t *s = session->md_rng;
  const uint64_t result =
      ((s[0] + s[3]) << 23 | (s[0] + s[3]) >> (64 - 23)) + s[0];
  const uint64_t t = s[1] << 17;
  s[2] ^= s[0];
  s[3] ^= s[1];
  s[1] ^= s[2];
  s[0] ^= s[3];
  s[2] ^= t;
  s[3] = (s[3] << 45) | (s[3] >> (64 - 45));
  return result;
}

double md_uniform(LammpcSession *session) {
  return (double)(md_random_u64(session) >> 11) *
         (1.0 / 9007199254740992.0);
}

double md_normal(LammpcSession *session) {
  const double u1 = std::max(md_uniform(session), 1.0e-16);
  const double u2 = md_uniform(session);
  return std::sqrt(-2.0 * std::log(u1)) *
         std::cos(6.283185307179586476925286766559 * u2);
}

double md_plane_denominator(const LammpcSession *session) {
  double denominator = 0.0;
  for (size_t i = 0; i < session->md_normal.size(); ++i) {
    const double normal = session->md_normal[i];
    denominator += normal * normal / session->md_masses[i];
  }
  return denominator;
}

void md_restore_position(LammpcSession *session) {
  double distance = 0.0;
  if (session->md_config.coordinates == LAMMPC_PLANE_MASS_WEIGHTED) {
    for (size_t i = 0; i < session->md_positions.size(); ++i)
      distance += session->md_normal[i] *
                  (std::sqrt(session->md_masses[i]) *
                       session->md_positions[i] -
                   session->md_point[i]);
    for (size_t i = 0; i < session->md_positions.size(); ++i)
      session->md_positions[i] -=
          session->md_normal[i] * distance /
          std::sqrt(session->md_masses[i]);
    return;
  }
  for (size_t i = 0; i < session->md_positions.size(); ++i)
    distance += session->md_normal[i] *
                (session->md_positions[i] - session->md_point[i]);
  const double coefficient = distance / md_plane_denominator(session);
  for (size_t i = 0; i < session->md_positions.size(); ++i)
    session->md_positions[i] -=
        session->md_normal[i] * coefficient / session->md_masses[i];
}

void md_project_velocity(LammpcSession *session) {
  std::fill(session->md_removed_velocities.begin(),
            session->md_removed_velocities.end(), 0.0);
  double component = 0.0;
  if (session->md_config.coordinates == LAMMPC_PLANE_MASS_WEIGHTED) {
    for (size_t i = 0; i < session->md_velocities.size(); ++i)
      component += session->md_normal[i] *
                   std::sqrt(session->md_masses[i]) *
                   session->md_velocities[i];
    for (size_t i = 0; i < session->md_velocities.size(); ++i) {
      const double removed =
          session->md_normal[i] * component /
          std::sqrt(session->md_masses[i]);
      session->md_velocities[i] -= removed;
      session->md_removed_velocities[i] = removed;
    }
    return;
  }
  for (size_t i = 0; i < session->md_velocities.size(); ++i)
    component += session->md_normal[i] * session->md_velocities[i];
  const double coefficient = component / md_plane_denominator(session);
  for (size_t i = 0; i < session->md_velocities.size(); ++i) {
    const double removed = session->md_normal[i] * coefficient /
                           session->md_masses[i];
    session->md_velocities[i] -= removed;
    session->md_removed_velocities[i] = removed;
  }
}

void md_project_force(LammpcSession *session) {
  if (session->md_config.coordinates == LAMMPC_PLANE_MASS_WEIGHTED) {
    double component = 0.0;
    for (size_t i = 0; i < session->md_forces.size(); ++i)
      component += session->md_normal[i] * session->md_forces[i] /
                   std::sqrt(session->md_masses[i]);
    session->md_constraint_force = component;
    for (size_t i = 0; i < session->md_forces.size(); ++i)
      session->md_forces[i] -=
          std::sqrt(session->md_masses[i]) * session->md_normal[i] *
          component;
    return;
  }
  double normal_force = 0.0;
  double acceleration = 0.0;
  for (size_t i = 0; i < session->md_forces.size(); ++i) {
    normal_force += session->md_normal[i] * session->md_forces[i];
    acceleration += session->md_normal[i] * session->md_forces[i] /
                    session->md_masses[i];
  }
  session->md_constraint_force = normal_force;
  const double coefficient = acceleration / md_plane_denominator(session);
  for (size_t i = 0; i < session->md_forces.size(); ++i)
    session->md_forces[i] -= session->md_normal[i] * coefficient;
}

int md_evaluate(LammpcSession *session) {
  std::vector<double> forces_h_bohr(session->md_positions.size(), 0.0);
  double energy_h = 0.0;
  double stress_h_bohr3[9] = {0};
  if (evaluate(session, session->n_atoms, session->md_positions.data(),
               &energy_h, forces_h_bohr.data(), stress_h_bohr3) != 0)
    return -1;
  for (size_t i = 0; i < forces_h_bohr.size(); ++i)
    session->md_forces[i] =
        forces_h_bohr[i] * kEvPerHartree / kAngstromPerBohr;
  ++session->md_force_evaluations;
  md_project_force(session);
  return 0;
}

int nve_coordinate_count(const LammpcSession *session, size_t *coordinate_count) {
  if (!session || !session->lmp || session->n_atoms == 0 ||
      session->atom_types.size() != session->n_atoms) {
    store_error("NVE state requires a live session geometry");
    return -1;
  }
  *coordinate_count = session->n_atoms * 3u;
  return 0;
}

int nve_sync_state(LammpcSession *session) {
  size_t coordinate_count = 0;
  if (nve_coordinate_count(session, &coordinate_count) != 0)
    return -1;
  session->nve_positions.resize(coordinate_count, 0.0);
  session->nve_velocities.resize(coordinate_count, 0.0);
  lammps_gather_atoms(session->lmp, "x", 1, 3, session->nve_positions.data());
  if (lammps_has_error(session->lmp)) {
    char buffer[512] = {0};
    if (lammps_get_last_error_message(session->lmp, buffer, (int)sizeof(buffer)) > 0)
      store_error(buffer);
    else
      store_error("lammps_gather_atoms(x) failed");
    return -1;
  }
  lammps_gather_atoms(session->lmp, "v", 1, 3, session->nve_velocities.data());
  if (lammps_has_error(session->lmp)) {
    char buffer[512] = {0};
    if (lammps_get_last_error_message(session->lmp, buffer, (int)sizeof(buffer)) > 0)
      store_error(buffer);
    else
      store_error("lammps_gather_atoms(v) failed");
    return -1;
  }
  session->nve_observation.potential_energy_ev =
      lammps_get_thermo(session->lmp, "pe");
  session->nve_observation.kinetic_energy_ev =
      lammps_get_thermo(session->lmp, "ke");
  return 0;
}

int md_apply_hard_boundary(LammpcSession *session,
                           const std::vector<double> &previous_positions,
                           const std::vector<double> &previous_velocities) {
  const HardBoundaryHook hook = session->hard_boundary;
  if (!hook.fn)
    return 0;

  auto *lmp = static_cast<LAMMPS_NS::LAMMPS *>(session->lmp);
  LammpcBoundaryBox box = {};
  for (int dim = 0; dim < 3; ++dim) {
    box.lo[dim] = lmp->domain->boxlo[dim];
    box.hi[dim] = lmp->domain->boxhi[dim];
    box.periodic[dim] = lmp->domain->periodicity[dim];
  }
  box.xy = lmp->domain->xy;
  box.xz = lmp->domain->xz;
  box.yz = lmp->domain->yz;
  box.triclinic = lmp->domain->triclinic;

  const int decision = lammpc_hard_boundary_apply(
      hook.fn, hook.user, static_cast<int>(session->n_atoms),
      session->md_positions.data(), session->md_velocities.data(),
      previous_positions.data(), previous_velocities.data(), &box);
  if (decision < 0) {
    store_error("hard boundary callback failed");
    return -1;
  }
  return decision;
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

void forget_hard_boundary(LammpcSession *session) {
  if (!session)
    return;
  session->hard_boundary = {};
  if (session->lmp)
    lammpc_embed_boundary_disarm(session->lmp);
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
  forget_hard_boundary(session);
  destroy_instance(session);
  delete session;
}

int lammpc_session_md_initialize(
    LammpcSession *session, const void *force_input_capnp,
    size_t force_input_capnp_size_bytes, const double *point,
    const double *normal, size_t coordinate_count,
    const LammpcMdConfig *config) {
  if (!session || !point || !normal || !config) {
    store_error("md_initialize received a null argument");
    return -1;
  }
  if (!(std::isfinite(config->time_step_ps) && config->time_step_ps > 0.0 &&
        std::isfinite(config->thermal_energy_ev) &&
        config->thermal_energy_ev > 0.0 &&
        std::isfinite(config->friction_per_ps) &&
        config->friction_per_ps > 0.0)) {
    store_error("constrained MD controls must be finite and positive");
    return -1;
  }
  if (config->coordinates != LAMMPC_PLANE_MASS_WEIGHTED &&
      config->coordinates != LAMMPC_PLANE_CARTESIAN) {
    store_error("unknown constrained MD coordinate convention");
    return -1;
  }

  struct capn arena;
  ForceInput_ptr force_input;
  if (lammpc_force_input_root(force_input_capnp,
                              force_input_capnp_size_bytes, &arena,
                              &force_input) != 0) {
    store_error("invalid ForceInput message");
    return -1;
  }
  size_t n_atoms = 0;
  int has_cell = 0;
  if (lammpc_force_input_atom_count(force_input, &n_atoms, &has_cell) != 0 ||
      coordinate_count != n_atoms * 3u) {
    lammpc_params_release(&arena);
    store_error("plane coordinate count differs from ForceInput geometry");
    return -1;
  }
  std::vector<double> positions(coordinate_count, 0.0);
  std::vector<int> atomic_numbers(n_atoms, 0);
  double cell_ang[9] = {0};
  if (lammpc_force_input_copy_geometry(
          force_input, positions.data(), atomic_numbers.data(), n_atoms,
          cell_ang, &has_cell) != 0) {
    lammpc_params_release(&arena);
    store_error("invalid ForceInput geometry");
    return -1;
  }
  lammpc_params_release(&arena);

  if (!geometry_matches(session, n_atoms, cell_ang, has_cell) &&
      build_instance(session, n_atoms, positions.data(), atomic_numbers.data(),
                     cell_ang, has_cell) != 0)
    return -1;

  session->md_config = *config;
  session->md_point.assign(point, point + coordinate_count);
  session->md_normal.assign(normal, normal + coordinate_count);
  session->md_positions = positions;
  session->md_velocities.assign(coordinate_count, 0.0);
  session->md_forces.assign(coordinate_count, 0.0);
  session->md_removed_velocities.assign(coordinate_count, 0.0);
  session->md_masses.assign(coordinate_count, 0.0);
  double normal_norm2 = 0.0;
  for (size_t i = 0; i < coordinate_count; ++i) {
    if (!std::isfinite(point[i]) || !std::isfinite(normal[i])) {
      store_error("plane point and normal must be finite");
      return -1;
    }
    normal_norm2 += normal[i] * normal[i];
    const size_t atom = i / 3u;
    const int type = session->atom_types[atom];
    const double mass =
        (type > 0 && (size_t)type <= session->params.masses.size())
            ? session->params.masses[(size_t)type - 1u]
            : 1.0;
    if (!(std::isfinite(mass) && mass > 0.0)) {
      store_error("constrained MD needs finite positive masses");
      return -1;
    }
    session->md_masses[i] = mass;
  }
  if (!(normal_norm2 > 0.0)) {
    store_error("plane normal vanishes");
    return -1;
  }
  const double inverse_normal_norm = 1.0 / std::sqrt(normal_norm2);
  for (double &value : session->md_normal)
    value *= inverse_normal_norm;
  if (config->coordinates == LAMMPC_PLANE_CARTESIAN &&
      !(md_plane_denominator(session) > 0.0)) {
    store_error("Cartesian plane has a zero mass-metric denominator");
    return -1;
  }

  md_seed(session, config->seed);
  for (size_t i = 0; i < coordinate_count; ++i) {
    const double sigma =
        std::sqrt(config->thermal_energy_ev /
                  (session->md_masses[i] * kMassVelocityToEv));
    session->md_velocities[i] = sigma * md_normal(session);
  }
  md_restore_position(session);
  md_project_velocity(session);
  session->md_force_evaluations = 0;
  session->md_steps = 0;
  session->md_ready = true;
  if (md_evaluate(session) != 0) {
    session->md_ready = false;
    return -1;
  }
  clear_error();
  return 0;
}

int lammpc_session_md_set_velocities(
    LammpcSession *session, const double *velocities_ang_ps,
    size_t coordinate_count, unsigned long long seed) {
  if (!session || !session->md_ready || !velocities_ang_ps ||
      coordinate_count != session->md_velocities.size()) {
    store_error("constrained MD velocity input has the wrong size");
    return -1;
  }
  for (size_t i = 0; i < coordinate_count; ++i) {
    if (!std::isfinite(velocities_ang_ps[i])) {
      store_error("constrained MD velocities must be finite");
      return -1;
    }
  }
  std::copy(velocities_ang_ps, velocities_ang_ps + coordinate_count,
            session->md_velocities.begin());
  md_project_velocity(session);
  md_seed(session, seed);
  clear_error();
  return 0;
}

int lammpc_session_masses(LammpcSession *session, double *masses_amu,
                          size_t atom_count) {
  if (!session || !masses_amu || atom_count != session->n_atoms ||
      session->atom_types.size() != atom_count || !session->lmp) {
    store_error("session mass output has the wrong size");
    return -1;
  }
  auto *atom_data = static_cast<LAMMPS_NS::LAMMPS *>(session->lmp)->atom;
  for (size_t atom = 0; atom < atom_count; ++atom) {
    const int type = session->atom_types[atom];
    if (type <= 0 || type > atom_data->ntypes) {
      store_error("session atom type has no mass");
      return -1;
    }
    const double mass = atom_data->mass[type];
    if (!(std::isfinite(mass) && mass > 0.0)) {
      store_error("session mass must be finite and positive");
      return -1;
    }
    masses_amu[atom] = mass;
  }
  clear_error();
  return 0;
}

int lammpc_session_md_step(LammpcSession *session, size_t step_count) {
  if (!session || !session->md_ready) {
    store_error("constrained MD is not initialized");
    return -1;
  }
  const double h = session->md_config.time_step_ps;
  const double decay =
      std::exp(-session->md_config.friction_per_ps * h);
  std::vector<double> previous_positions;
  std::vector<double> previous_velocities;
  std::vector<double> previous_forces;
  std::vector<double> previous_removed_velocities;
  for (size_t step = 0; step < step_count; ++step) {
    const bool hard_boundary_armed = session->hard_boundary.fn != nullptr;
    previous_positions = session->md_positions;
    previous_velocities = session->md_velocities;
    previous_forces = session->md_forces;
    previous_removed_velocities = session->md_removed_velocities;
    const double previous_constraint_force = session->md_constraint_force;
    uint64_t previous_rng[4] = {};
    std::memcpy(previous_rng, session->md_rng, sizeof(previous_rng));
    auto rollback_step = [&]() {
      session->md_positions = previous_positions;
      session->md_velocities = previous_velocities;
      session->md_forces = previous_forces;
      session->md_removed_velocities = previous_removed_velocities;
      session->md_constraint_force = previous_constraint_force;
      std::memcpy(session->md_rng, previous_rng, sizeof(previous_rng));
      if (!session->lmp ||
          static_cast<LAMMPS_NS::LAMMPS *>(session->lmp)->atom->nlocal !=
              static_cast<int>(session->n_atoms)) {
        session->md_ready = false;
      }
    };
    for (size_t i = 0; i < session->md_velocities.size(); ++i)
      session->md_velocities[i] +=
          0.5 * h * session->md_forces[i] /
          (session->md_masses[i] * kMassVelocityToEv);
    for (size_t i = 0; i < session->md_positions.size(); ++i)
      session->md_positions[i] += 0.5 * h * session->md_velocities[i];
    md_restore_position(session);
    for (size_t i = 0; i < session->md_velocities.size(); ++i) {
      const double kick =
          std::sqrt(session->md_config.thermal_energy_ev *
                    (1.0 - decay * decay) /
                    (session->md_masses[i] * kMassVelocityToEv));
      session->md_velocities[i] =
          decay * session->md_velocities[i] + kick * md_normal(session);
    }
    md_project_velocity(session);
    for (size_t i = 0; i < session->md_positions.size(); ++i)
      session->md_positions[i] += 0.5 * h * session->md_velocities[i];
    md_restore_position(session);
    if (md_evaluate(session) != 0) {
      rollback_step();
      return -1;
    }
    for (size_t i = 0; i < session->md_velocities.size(); ++i)
      session->md_velocities[i] +=
          0.5 * h * session->md_forces[i] /
          (session->md_masses[i] * kMassVelocityToEv);
    md_project_velocity(session);
    if (hard_boundary_armed) {
      const int decision = md_apply_hard_boundary(
          session, previous_positions, previous_velocities);
      if (decision < 0) {
        rollback_step();
        return -1;
      }
      if (decision == LAMMPC_BOUNDARY_REJECT) {
        session->md_forces = previous_forces;
        session->md_removed_velocities = previous_removed_velocities;
        session->md_constraint_force = previous_constraint_force;
      } else if (decision == LAMMPC_BOUNDARY_REFLECT) {
        md_restore_position(session);
        md_project_velocity(session);
        if (md_evaluate(session) != 0) {
          rollback_step();
          return -1;
        }
      }
    }
    ++session->md_steps;
  }
  clear_error();
  return 0;
}

int lammpc_session_md_state(
    LammpcSession *session, double *positions_ang,
    double *velocities_ang_ps, double *projected_forces_ev_ang,
    double *removed_velocities_ang_ps, size_t coordinate_count,
    LammpcMdObservation *observation) {
  if (!session || !session->md_ready ||
      coordinate_count != session->md_positions.size()) {
    store_error("constrained MD state request has the wrong size");
    return -1;
  }
  if (positions_ang)
    std::copy(session->md_positions.begin(), session->md_positions.end(),
              positions_ang);
  if (velocities_ang_ps)
    std::copy(session->md_velocities.begin(), session->md_velocities.end(),
              velocities_ang_ps);
  if (projected_forces_ev_ang)
    std::copy(session->md_forces.begin(), session->md_forces.end(),
              projected_forces_ev_ang);
  if (removed_velocities_ang_ps)
    std::copy(session->md_removed_velocities.begin(),
              session->md_removed_velocities.end(),
              removed_velocities_ang_ps);
  if (observation) {
    observation->constraint_force = session->md_constraint_force;
    observation->force_evaluations = session->md_force_evaluations;
    observation->steps = session->md_steps;
  }
  clear_error();
  return 0;
}

int lammpc_session_nve_set_positions(LammpcSession *session,
                                     const double *positions_ang,
                                     size_t coordinate_count,
                                     const double *cell_ang, size_t cell_len) {
  if (!session || !positions_ang) {
    store_error("nve_set_positions received a null argument");
    return -1;
  }
  size_t expected_coordinates = 0;
  if (nve_coordinate_count(session, &expected_coordinates) != 0)
    return -1;
  if (coordinate_count != expected_coordinates) {
    store_error("nve_set_positions coordinate count mismatch");
    return -1;
  }
  if ((cell_ang && cell_len != 9u) || (!cell_ang && cell_len != 0u)) {
    store_error("nve_set_positions box must have exactly 9 entries");
    return -1;
  }
  if (cell_ang) {
    if (std::fabs(cell_ang[1]) > 1e-12 || std::fabs(cell_ang[2]) > 1e-12 ||
        std::fabs(cell_ang[5]) > 1e-12) {
      store_error("nve_set_positions supports restricted triclinic boxes only");
      return -1;
    }
    if (!(cell_ang[0] > 0.0 && cell_ang[4] > 0.0 && cell_ang[8] > 0.0)) {
      store_error("nve_set_positions box diagonal entries must be positive");
      return -1;
    }
    const std::string cmd =
        "change_box all x final 0 " + format_double17(cell_ang[0]) +
        " y final 0 " + format_double17(cell_ang[4]) + " z final 0 " +
        format_double17(cell_ang[8]) + " xy final " + format_double17(cell_ang[3]) +
        " xz final " + format_double17(cell_ang[6]) + " yz final " +
        format_double17(cell_ang[7]) + " remap units box";
    if (run_cmd(session->lmp, cmd) != 0)
      return -1;
    session->has_cell = 1;
    std::memcpy(session->cell_ang, cell_ang, sizeof(session->cell_ang));
  }
  lammps_scatter_atoms(session->lmp, "x", 1, 3, (void *)positions_ang);
  if (lammps_has_error(session->lmp)) {
    char buffer[512] = {0};
    if (lammps_get_last_error_message(session->lmp, buffer, (int)sizeof(buffer)) > 0)
      store_error(buffer);
    else
      store_error("lammps_scatter_atoms(x) failed");
    return -1;
  }
  if (run_cmd(session->lmp, "run 0 pre yes post no") != 0)
    return -1;
  if (nve_sync_state(session) != 0)
    return -1;
  clear_error();
  return 0;
}

int lammpc_session_nve_set_velocities(LammpcSession *session,
                                      const double *velocities_ang_ps,
                                      size_t coordinate_count) {
  if (!session || !velocities_ang_ps) {
    store_error("nve_set_velocities received a null argument");
    return -1;
  }
  size_t expected_coordinates = 0;
  if (nve_coordinate_count(session, &expected_coordinates) != 0)
    return -1;
  if (coordinate_count != expected_coordinates) {
    store_error("nve_set_velocities coordinate count mismatch");
    return -1;
  }
  lammps_scatter_atoms(session->lmp, "v", 1, 3, (void *)velocities_ang_ps);
  if (lammps_has_error(session->lmp)) {
    char buffer[512] = {0};
    if (lammps_get_last_error_message(session->lmp, buffer, (int)sizeof(buffer)) > 0)
      store_error(buffer);
    else
      store_error("lammps_scatter_atoms(v) failed");
    return -1;
  }
  if (run_cmd(session->lmp, "run 0 pre yes post no") != 0)
    return -1;
  if (nve_sync_state(session) != 0)
    return -1;
  clear_error();
  return 0;
}

int lammpc_session_nve_step(LammpcSession *session, size_t step_count,
                            double time_step_ps) {
  if (!session) {
    store_error("nve_step received a null session");
    return -1;
  }
  if (!(std::isfinite(time_step_ps) && time_step_ps > 0.0)) {
    store_error("nve_step needs a finite positive timestep");
    return -1;
  }
  size_t expected_coordinates = 0;
  if (nve_coordinate_count(session, &expected_coordinates) != 0)
    return -1;
  (void)expected_coordinates;
  if (!session->nve_fix_installed) {
    if (run_cmd(session->lmp, "fix __lammpc_nve all nve") != 0)
      return -1;
    session->nve_fix_installed = true;
  }
  if (run_cmd(session->lmp, "timestep " + format_double17(time_step_ps)) != 0)
    return -1;
  if (step_count > 0 &&
      run_cmd(session->lmp, "run " + std::to_string(step_count) +
                                 " pre yes post no") != 0)
    return -1;
  if (step_count == 0 && run_cmd(session->lmp, "run 0 pre yes post no") != 0)
    return -1;
  if (nve_sync_state(session) != 0)
    return -1;
  clear_error();
  return 0;
}

int lammpc_session_nve_state(LammpcSession *session, double *positions_ang,
                             double *velocities_ang_ps,
                             size_t coordinate_count,
                             LammpcNveObservation *observation) {
  if (!session) {
    store_error("nve_state received a null session");
    return -1;
  }
  size_t expected_coordinates = 0;
  if (nve_coordinate_count(session, &expected_coordinates) != 0)
    return -1;
  if (coordinate_count != expected_coordinates) {
    store_error("nve_state coordinate count mismatch");
    return -1;
  }
  if (nve_sync_state(session) != 0)
    return -1;
  if (positions_ang) {
    std::copy(session->nve_positions.begin(), session->nve_positions.end(),
              positions_ang);
  }
  if (velocities_ang_ps) {
    std::copy(session->nve_velocities.begin(), session->nve_velocities.end(),
              velocities_ang_ps);
  }
  if (observation)
    *observation = session->nve_observation;
  clear_error();
  return 0;
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

int lammpc_session_set_hard_boundary(LammpcSession *session,
                                     LammpcHardBoundaryFn fn, void *user) {
  if (!session) {
    store_error("hard boundary requires a session");
    return -1;
  }
  if (!fn) {
    session->hard_boundary = {};
    if (session->lmp && lammpc_embed_boundary_detach(session->lmp) != 0) {
      store_error(lammpc_embed_boundary_error());
      return -1;
    }
    clear_error();
    return 0;
  }
  session->hard_boundary = HardBoundaryHook{fn, user};
  clear_error();
  return 0;
}

int lammpc_session_clear_hard_boundary(LammpcSession *session) {
  return lammpc_session_set_hard_boundary(session, nullptr, nullptr);
}

int lammpc_session_integrate(LammpcSession *session, int nsteps) {
  if (!session || nsteps < 0) {
    store_error("integrate needs a session and a non-negative step count");
    return -1;
  }
  if (!session->lmp) {
    store_error("integrate needs a live LAMMPS instance");
    return -1;
  }

  /* Unset: run the steps the host already configured. No fix, no snapshot. */
  const HardBoundaryHook hook = session->hard_boundary;
  const int armed = hook.fn != nullptr;
  if (armed &&
      lammpc_embed_boundary_attach(session->lmp, hook.fn, hook.user) != 0) {
    store_error(lammpc_embed_boundary_error());
    return -1;
  }

  char command[64];
  std::snprintf(command, sizeof(command), "run %d", nsteps);
  const int ran = run_cmd(session->lmp, command);
  if (armed) {
    const int detached = lammpc_embed_boundary_detach(session->lmp);
    if (ran == 0 && detached != 0) {
      store_error(lammpc_embed_boundary_error());
      return -1;
    }
  }
  if (ran != 0)
    return -1;
  clear_error();
  return 0;
}

} // extern "C"
