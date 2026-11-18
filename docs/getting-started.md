# Getting started

This page builds the stub library and runs a program that inspects its public
ABI. The stub exports the complete ABI but does not link LAMMPS.

## Requirements

Install Pixi. Pixi creates the C and C++ build environments declared in
`pixi.toml`.

## Build the stub

From the repository root, run:

```sh
pixi run test-stub
```

This command configures Meson in `build`, compiles `liblammpc.so`, and runs
the C test suite. The installed header set consists of `lammpc.h`,
`lammpc_boundary.h`, and `lammpc_features.h`.

The `public-api-example` test builds and runs `examples/inspect.c`. Run it
directly after the build with:

```sh
LD_LIBRARY_PATH=build build/public_api_example
```

The program prints the library version, ABI generation, runtime availability,
and number of feature entries. It also checks the documented size query for a
`Capabilities` message.

## Build with LAMMPS

The embed environment downloads and compiles the pinned LAMMPS release:

```sh
pixi run -e embed build-lammps
pixi run -e embed test-embed
```

The install is stored in `.lammps-install`. The second command configures
`build-embed` with `with_lammps=true`, compiles lammpc, and runs evaluations
against the real library.

The LAMMPS configuration is serial. It disables MPI, OpenMP, and build tools.
It enables the MANYBODY package for the EAM NVE test.

## Check the documentation

Run:

```sh
pixi run docs
```

This command builds Doxygen with warnings treated as errors and checks every
external documentation link. DOI links are checked through the DOI handle
API.

## Select the correct library

Call `lammpc_available()` after loading the shared library. A return value of
one identifies the LAMMPS embed. A return value of zero identifies the stub.
Both variants expose the same symbols. Evaluation calls on the stub fail and
explain the failure in `LammpcResult.message` or `lammpc_last_error()`.
