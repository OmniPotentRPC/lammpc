/**
 * fix lammpc/boundary, registered on the embed's LAMMPS instance.
 *
 * This is not an upstream LAMMPS style. The creator is inserted into that
 * instance's fix map, and the fix is added only for a run whose session
 * has a host callback. end_of_step() sees positions, velocities, and the
 * box after the integrator (and LAMMPS PBC) have finished the step.
 * Reject writes the previous positions, velocities, and image flags back.
 * Reflect writes the host buffers back. Accept leaves the integrator state.
 */
#include "lammpc_boundary_embed.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "lammps.h"
#include "modify.h"

#include <cstdio>
#include <cstring>
#include <exception>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

using LAMMPS_NS::bigint;
using LAMMPS_NS::imageint;
using LAMMPS_NS::tagint;

thread_local char g_boundary_error[512] = "";

struct BoundaryBinding {
  LammpcHardBoundaryFn fn = nullptr;
  void *user = nullptr;
};

std::unordered_map<LAMMPS_NS::LAMMPS *, BoundaryBinding> g_bindings;

void set_boundary_error(const char *message) {
  if (!message)
    message = "lammpc hard boundary failed";
  std::snprintf(g_boundary_error, sizeof(g_boundary_error), "%s", message);
}

class FixLammpcBoundary : public LAMMPS_NS::Fix {
 public:
  FixLammpcBoundary(LAMMPS_NS::LAMMPS *lmp, int narg, char **arg)
      : Fix(lmp, narg, arg) {
    nevery = 1;
    time_integrate = 0;
    dynamic_group_allow = 1;
  }

  int setmask() override { return LAMMPS_NS::FixConst::END_OF_STEP; }

  void setup(int) override { capture_previous(); }

  void end_of_step() override {
    auto found = g_bindings.find(lmp);
    if (found == g_bindings.end() || !found->second.fn)
      return;

    const int n_atoms = serial_atom_count();
    if (previous_images_.size() != static_cast<size_t>(n_atoms) ||
        previous_positions_.size() != static_cast<size_t>(n_atoms) * 3u)
      error->all(FLERR, "lammpc hard boundary has no previous state");
    std::vector<double> positions(static_cast<size_t>(n_atoms) * 3u);
    std::vector<double> velocities(static_cast<size_t>(n_atoms) * 3u);
    std::vector<int> local_of_tag(static_cast<size_t>(n_atoms), -1);
    read_atoms(n_atoms, local_of_tag, positions.data(), velocities.data());

    LammpcBoundaryBox box = {};
    fill_box(&box);
    const int decision = lammpc_hard_boundary_apply(
        found->second.fn, found->second.user, n_atoms, positions.data(),
        velocities.data(), previous_positions_.data(),
        previous_velocities_.data(), &box);
    if (decision < 0) {
      write_atoms(n_atoms, local_of_tag, positions.data(), velocities.data(),
                  1);
      capture_previous();
      error->all(FLERR, "lammpc hard boundary callback failed");
    }
    if (decision == LAMMPC_BOUNDARY_REJECT)
      write_atoms(n_atoms, local_of_tag, positions.data(), velocities.data(),
                  1);
    else if (decision == LAMMPC_BOUNDARY_REFLECT)
      write_atoms(n_atoms, local_of_tag, positions.data(), velocities.data(),
                  0);
    capture_previous();
  }

 private:
  std::vector<double> previous_positions_;
  std::vector<double> previous_velocities_;
  std::vector<imageint> previous_images_;

  int serial_atom_count() const {
    if (static_cast<bigint>(atom->nlocal) != atom->natoms)
      error->all(FLERR,
                 "lammpc hard boundary requires a serial atom decomposition");
    if (atom->natoms > static_cast<bigint>(0x7fffffff))
      error->all(FLERR, "lammpc hard boundary atom count exceeds int");
    return atom->nlocal;
  }

  void fill_box(LammpcBoundaryBox *box) const {
    for (int dim = 0; dim < 3; ++dim) {
      box->lo[dim] = domain->boxlo[dim];
      box->hi[dim] = domain->boxhi[dim];
      box->periodic[dim] = domain->periodicity[dim];
    }
    box->xy = domain->xy;
    box->xz = domain->xz;
    box->yz = domain->yz;
    box->triclinic = domain->triclinic;
  }

  void map_tags(int n_atoms, std::vector<int> *local_of_tag) const {
    local_of_tag->assign(static_cast<size_t>(n_atoms), -1);
    for (int local = 0; local < n_atoms; ++local) {
      const tagint tag = atom->tag[local];
      if (tag < 1 || tag > static_cast<tagint>(n_atoms))
        error->all(FLERR, "lammpc hard boundary requires tags 1..N");
      (*local_of_tag)[static_cast<size_t>(tag - 1)] = local;
    }
    for (int caller = 0; caller < n_atoms; ++caller) {
      if ((*local_of_tag)[static_cast<size_t>(caller)] < 0)
        error->all(FLERR, "lammpc hard boundary is missing an atom tag");
    }
  }

  void read_atoms(int n_atoms, std::vector<int> &local_of_tag, double *positions,
                  double *velocities) const {
    map_tags(n_atoms, &local_of_tag);
    for (int caller = 0; caller < n_atoms; ++caller) {
      const int local = local_of_tag[static_cast<size_t>(caller)];
      std::memcpy(positions + 3 * caller, atom->x[local], 3 * sizeof(double));
      std::memcpy(velocities + 3 * caller, atom->v[local], 3 * sizeof(double));
    }
  }

  void write_atoms(int n_atoms, const std::vector<int> &local_of_tag,
                   const double *positions, const double *velocities,
                   int restore_images) const {
    for (int caller = 0; caller < n_atoms; ++caller) {
      const int local = local_of_tag[static_cast<size_t>(caller)];
      std::memcpy(atom->x[local], positions + 3 * caller, 3 * sizeof(double));
      std::memcpy(atom->v[local], velocities + 3 * caller, 3 * sizeof(double));
      if (restore_images)
        atom->image[local] = previous_images_[static_cast<size_t>(caller)];
    }
  }

  void capture_previous() {
    const int n_atoms = serial_atom_count();
    previous_positions_.assign(static_cast<size_t>(n_atoms) * 3u, 0.0);
    previous_velocities_.assign(static_cast<size_t>(n_atoms) * 3u, 0.0);
    previous_images_.assign(static_cast<size_t>(n_atoms), 0);
    std::vector<int> local_of_tag;
    read_atoms(n_atoms, local_of_tag, previous_positions_.data(),
               previous_velocities_.data());
    for (int caller = 0; caller < n_atoms; ++caller) {
      const int local = local_of_tag[static_cast<size_t>(caller)];
      previous_images_[static_cast<size_t>(caller)] = atom->image[local];
    }
  }
};

LAMMPS_NS::Fix *create_boundary_fix(LAMMPS_NS::LAMMPS *lmp, int narg,
                                    char **arg) {
  return new FixLammpcBoundary(lmp, narg, arg);
}

}  // namespace

extern "C" {

int lammpc_embed_boundary_attach(void *handle, LammpcHardBoundaryFn fn,
                                 void *user) {
  if (!handle || !fn) {
    set_boundary_error("hard boundary attach needs an instance and a callback");
    return -1;
  }
  auto *lmp = static_cast<LAMMPS_NS::LAMMPS *>(handle);
  try {
    (*lmp->modify->fix_map)["lammpc/boundary"] = &create_boundary_fix;
    g_bindings[lmp] = BoundaryBinding{fn, user};
    if (lmp->modify->get_fix_by_id("lammpc_boundary") == nullptr)
      lmp->modify->add_fix("lammpc_boundary all lammpc/boundary", 0);
  } catch (const std::exception &ex) {
    g_bindings.erase(lmp);
    set_boundary_error(ex.what());
    return -1;
  } catch (...) {
    g_bindings.erase(lmp);
    set_boundary_error("hard boundary fix registration failed");
    return -1;
  }
  g_boundary_error[0] = '\0';
  return 0;
}

int lammpc_embed_boundary_detach(void *handle) {
  if (!handle)
    return 0;
  auto *lmp = static_cast<LAMMPS_NS::LAMMPS *>(handle);
  g_bindings.erase(lmp);
  try {
    if (lmp->modify->get_fix_by_id("lammpc_boundary") != nullptr)
      lmp->modify->delete_fix("lammpc_boundary");
  } catch (const std::exception &ex) {
    set_boundary_error(ex.what());
    return -1;
  } catch (...) {
    set_boundary_error("hard boundary fix removal failed");
    return -1;
  }
  g_boundary_error[0] = '\0';
  return 0;
}

void lammpc_embed_boundary_disarm(void *handle) {
  if (!handle)
    return;
  g_bindings.erase(static_cast<LAMMPS_NS::LAMMPS *>(handle));
}

const char *lammpc_embed_boundary_error(void) { return g_boundary_error; }

}  // extern "C"
