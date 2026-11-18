# lammpc

lammpc is a C ABI for evaluating LAMMPS interaction models in process. It
uses flat Cap'n Proto messages from
[potentials-schema](https://github.com/OmniPotentRPC/potentials-schema).
The library provides persistent sessions, projected Langevin dynamics,
unconstrained NVE dynamics, and host callbacks for hard boundaries.

lammpc currently supports LAMMPS `metal` units only. A `LammpsParams` message
that requests `real` units is rejected. See the
[unit reference](docs/reference.md#units) before preparing inputs.

## Documentation

The documentation follows four task and information types.

| Need | Page |
| --- | --- |
| Build and run the first program | [Getting started](docs/getting-started.md) |
| Evaluate forces or run dynamics | [How to use sessions](docs/how-to.md) |
| Understand the design and algorithms | [Design explanation](docs/explanation.md) |
| Look up units, errors, and functions | [C API reference](docs/reference.md) |

The generated Doxygen reference is built from `include/lammpc.h`,
`include/lammpc_boundary.h`, and `include/lammpc_features.h`.

## Build and test

Install Pixi, then run:

```sh
pixi run test-stub
pixi run test-stub-sanitize
pixi run docs
```

The real LAMMPS checks use the same commands as CI:

```sh
pixi run -e embed build-lammps
pixi run -e embed test-embed
pixi run -e embed test-embed-tsan
```

`build-lammps` installs the pinned serial LAMMPS release under
`.lammps-install`. The embed build links that installation and
[minimage](https://github.com/lode-org/minimage).

## Project layout

| Part | Path |
| --- | --- |
| Public C ABI | `include/` |
| Runnable public API example | `examples/inspect.c` |
| User documentation | `docs/` |
| Embed implementation | `src/lammpc_embed_api.cpp` |
| Stub implementation | `src/lammpc_stub.c` |
| Wire schema | `schema/Potentials.capnp` |
| Tests | `tests/` |
