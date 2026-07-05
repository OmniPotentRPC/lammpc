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

- One persistent LAMMPS instance per session (`lammpc_session_create_from_config`);
  the interaction model (units / pair_style / pair_coeff / type map) is
  parsed once from `LammpsParams`, never per step.
- Per-step geometry moves through native array scatter/gather
  (`lammps_scatter_atoms` / `lammps_gather_atoms` family), not command text.
- The step evaluation runs `run 0 pre yes post no`-equivalent minimal
  re-neighboring; no file traffic, no log parsing, no subprocess.
- The benchmark comparing both paths on the same pair style and geometry
  sequence is the merge gate for the embed shell.

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

The embed shell links liblammps (`-Dwith_lammps=true -Dlammps_root=...`)
through the LAMMPS C library interface (`library.h`).
