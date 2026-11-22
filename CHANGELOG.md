# Changelog

All notable changes are recorded here before a release.

## 0.2.0 (2026-11-22)

### Added

* Constrained-MD sessions can now replace and project an initial velocity state
  while restarting the deterministic Langevin stream, and expose the per-atom
  masses selected by the configured LAMMPS type map.
* Constrained MD accepts zero friction and integrates that choice as projected
  velocity Verlet. A zero-step call leaves an initialized session unchanged.

### Fixed

* Feature discovery lists lammpc_session_md_set_velocities and
  lammpc_session_masses.
* The hard-boundary callback now runs once after each completed constrained-MD
  step; rejected steps restore force bookkeeping and reflected steps are
  re-projected and force-re-evaluated, with failed step evaluations rolling back
  the pre-step state and RNG transactionally. Rejected Langevin proposals retain
  their consumed random draws so retries use fresh noise. Reflected, rejected,
  and invalid-callback LAMMPS integration steps all refresh forces, energy, and
  virial data at the resulting configuration. Callback state is owned by each
  session and fix instance, so independent sessions can set, clear, and destroy
  callbacks concurrently.

## 0.1.0 (2026-11-20)

### Added

* A stable C ABI with stub and embedded LAMMPS implementations.
* Persistent sessions for Cap'n Proto configuration and force evaluation.
* Projected BAOAB dynamics with Cartesian and mass weighted planes.
* Unconstrained NVE session dynamics through LAMMPS `fix nve`.
* Host hard boundary callbacks for direct integration and projected dynamics.
* Runtime capability and feature discovery.

### Fixed

* Periodic positions several cells away are folded with the full restricted
  triclinic lattice during session creation and NVE position updates.
* Hard boundary rejection restores complete accepted state while retaining
  consumed Langevin random draws.
* Hard boundary reflection refreshes forces at the reflected position.
* Independent session callbacks no longer share callback state.

### Documentation

* Added getting started, task, design, and C API reference pages.
* Added warning clean Doxygen and link checks to CI.
