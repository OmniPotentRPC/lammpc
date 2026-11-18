# lammpc

Stable C ABI over LAMMPS as a dynlib. rgpot (and through it eOn) drives
LAMMPS with direct in-process calls: flat binary Cap'n Proto messages both
ways, prefix-parameterized symbols (`lammpc_*`), no input scripts, no string
parsing in the force loop, no JSON. The same design as
[nwchemc](https://github.com/OmniPotentRPC/nwchemc) and
[cpmdc](https://github.com/OmniPotentRPC/cpmdc); the canonical wire contract
is [potentials-schema](https://github.com/OmniPotentRPC/potentials-schema)
(`LammpsParams` arm, `ForceInput` -> `PotentialResult`).

## Performance contract

The acceptance bar: the rgpot -> lammpc session path is FASTER per step than
eOn's current LAMMPS call style. The design decisions that carry that bar:

- One persistent LAMMPS instance per session
  (`lammpc_session_create_from_config`); the interaction model
  (units / pair style / coefficients / type map) is applied once from
  `LammpsParams` by calling the actual LAMMPS functions on the instance.
  Remaining setup text commands (region bounds, masses) use `%.17g` formatting.
- First-call geometry creation is full-precision and batched: periodic
  coordinates are wrapped through
  [minimage](https://github.com/lode-org/minimage), then all atoms are created
  in one `lammps_create_atoms` call. The build fails loudly if LAMMPS reports
  any atom-count mismatch (`requested != created || requested != present`).
- Per-step geometry writes straight into the instance's position arrays;
  forces read straight from the force arrays and energy from the pair
  accumulators after the direct pair compute call. No file traffic, no log
  parsing, no subprocess.
- The benchmark comparing both paths on the same pair style and geometry
  sequence (eOn checkout: `TheochemUI/eOn`) is the merge gate for the embed
  build.

## Layout

| Piece | Path |
| --- | --- |
| Public ABI | `include/lammpc.h`, `include/lammpc_features.h` |
| Stub build (no LAMMPS) | `src/lammpc_stub.c` |
| Cap'n Proto helpers | `src/lammpc_params.c` |
| Capability discovery | `src/lammpc_capabilities.c` |
| Vendored contract | `schema/Potentials.capnp` (byte-identical to the wrap pin; CI-enforced) |

## Minimum potential ABI profile

lammpc exports the full profile from potentials-schema `PROFILE.md`:
lifecycle (`lammpc_abi_version`, `lammpc_version`, `lammpc_last_error`,
`lammpc_available`, `lammpc_finalize`), configuration (`lammpc_configure`,
`lammpc_session_create_from_config`, `lammpc_session_configure`,
`lammpc_session_destroy`), evaluation (`lammpc_session_calculate_result`,
`lammpc_calculate_result_from_config`,
`lammpc_potential_result_size_for_force_input`), and discovery
(`lammpc_capabilities_result`, `lammpc_feature_*`). rgpot's profile-driven
loader (`rgpot::abi::ProfileLoader`) resolves the whole set from the prefix
`lammpc` with zero backend-specific loader code.

`CommonMethodSpec` overlay: a classical force field lowers none of the
electronic-structure fields, so `Capabilities.loweredCommonFields` is empty
and `lammpc_configure` rejects any set overlay field through
`lammpc_last_error` instead of silently ignoring it.

## Build

```sh
pixi run test-stub          # stub + cmocka suites
pixi run test-stub-sanitize # same under ASan/UBSan
```

The embed build links liblammps (`-Dwith_lammps=true -Dlammps_root=...`) and
calls its actual functions directly on the force loop. It also links minimage
for periodic wrap/face handling on first-call atom creation.

## Constrained molecular dynamics

The additive `lammpc_session_md_*` API runs projected BAOAB dynamics on a
configured session. `lammpc_session_md_initialize` receives one `ForceInput`,
a hyperplane, and metal-unit time step, thermal energy, friction, and seed.
`lammpc_session_md_step` advances an exact number of steps.
`lammpc_session_md_state` returns Cartesian positions, projected forces,
removed velocities, the scalar constraint force, and step and force-call
counts.

Both mass-weighted and Cartesian planes are supported. Position is restored
after each drift, force and velocity are projected every step, and Langevin
noise is projected onto the plane. The interface is independent of PAFI,
OH-TST, pair styles, and material types.

## Unconstrained NVE session dynamics

An additive `lammpc_session_nve_*` ABI provides an opt-in unconstrained path
for trajectory shooting and reversible dynamics checks:

- `lammpc_session_nve_set_positions`: fold far periodic images along periodic
  cell axes (leaving non-periodic axes untouched), scatter Cartesian positions
  on the live geometry, and optionally update the box.
- `lammpc_session_nve_set_velocities`: scatter Cartesian velocities.
- `lammpc_session_nve_step`: run exact velocity-Verlet steps with LAMMPS's own
  `fix nve` integrator (`run N`, no projection/restoration).
- `lammpc_session_nve_state`: gather Cartesian positions, velocities, and
  potential/kinetic energy from thermo.

This path does not modify the constrained BAOAB interface; both APIs coexist
on the same session type.
