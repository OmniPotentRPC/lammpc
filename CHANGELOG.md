# Changelog

All notable changes are recorded here before a release.

## Unreleased

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
