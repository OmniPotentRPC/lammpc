# C API reference

This page summarizes the installed headers. Doxygen provides the declaration
level reference. lammpc has no Rust API, Python binding, or command line
interface.

## Units

### LAMMPS configuration

| `LammpsParams.unitsStyle` | Result |
| --- | --- |
| Empty | Uses `metal` |
| `metal` | Accepted |
| `real` | Rejected with `only unitsStyle=metal is supported` |
| Any other value | Rejected with the same error |

`atomStyle` defaults to `atomic`. No other atom style is accepted.

### Evaluation messages

`ForceInput.lengthUnit` accepts angstrom, bohr, nanometer, and the aliases
implemented in `src/lammpc_params.c`. It defaults to angstrom.
`ForceInput.energyUnit` accepts Hartree, electronvolt, Rydberg, kilojoules per
mole, and kilocalories per mole. It defaults to electronvolt in the schema.
The returned `PotentialResult` uses the requested energy unit. Forces use that
energy unit divided by the requested length unit. Stress uses that energy unit
divided by the cube of the requested length unit.

`LammpcResult.energy_h` always uses Hartree. The legacy raw array functions
name their units directly: positions use angstroms, energy uses Hartree, and
forces or gradients use Hartree per Bohr. Those raw array calls are retained
for ABI discovery but the current embed requires session messages and returns
an error.

### Dynamics

| Quantity | Unit |
| --- | --- |
| Position and cell | angstrom |
| Velocity | angstrom per picosecond |
| Projected force | electronvolt per angstrom |
| NVE energy | electronvolt |
| MD time step | picosecond |
| MD friction | reciprocal picosecond |
| MD thermal energy | electronvolt |
| Mass | atomic mass unit |

## Return conventions

Most functions returning `int` return zero on success and minus one on
failure. The hard-boundary helper `lammpc_hard_boundary_apply()` returns a
decision code on success: accept, reject, or reflect.

Read the thread local `lammpc_last_error()` immediately after a failure from an
`int`-returning API. Successful embed calls clear that string. The stub keeps
a fixed unavailable message.

Evaluation functions return `LammpcResult`. `ok` is one on success and zero
on failure. `message` is always null terminated. On success, it contains
`ok`. On failure, it contains the relevant diagnostic.

Creation functions return `NULL` on failure. The result size query returns
zero for an invalid `ForceInput` or an unrepresentable size.

No C function throws a C++ exception across the ABI. As with ordinary C,
passing a dangling pointer or a buffer smaller than its declared size is
outside the contract.

## Lifecycle and discovery

| Function | Contract |
| --- | --- |
| `lammpc_version` | Returns library owned static text. |
| `lammpc_abi_version` | Returns `LAMMPC_ABI_VERSION`, currently zero. |
| `lammpc_available` | Returns one for the embed and zero for the stub. |
| `lammpc_last_error` | Returns thread local diagnostic text. |
| `lammpc_finalize` | No operation in the current ABI. |
| `lammpc_capabilities_result` | Writes the backend `Capabilities` message. |
| `lammpc_feature_count` | Returns the number of static feature entries. |
| `lammpc_feature_table` | Returns a library owned static table. |
| `lammpc_feature_find` | Returns a matching static entry or `NULL`. |

For a capabilities size query, pass `NULL`, zero, and a valid size output
pointer. The function returns minus one and writes the required byte count.
Repeat the call with a buffer of at least that size. A missing size output
pointer is an error.

Feature table pointers and their strings remain valid for the life of the
loaded library. Callers must not free or modify them.

## Configuration and session lifetime

| Function | Contract |
| --- | --- |
| `lammpc_set_params` | Legacy global entry point. It is unsupported by the embed because global parameters are disabled. |
| `lammpc_configure` | Validates a `PotentialConfig` by creating and destroying a temporary session. |
| `lammpc_session_create` | Creates an idle session from `LammpsParams`. |
| `lammpc_session_create_from_config` | Creates an idle session from the `lammps` arm. |
| `lammpc_session_set_params` | Replaces parameters and destroys live geometry. |
| `lammpc_session_configure` | Replaces configuration and destroys live geometry. |
| `lammpc_session_destroy` | Releases callback, LAMMPS, and session state. `NULL` is accepted. |

The parameter defaults are `metal`, `atomic`, and `p p p`. `pairStyle` and
`typeToAtomicNumber` are required. `masses` defaults each omitted type to one
atomic mass unit in the current implementation. Supply explicit positive
masses for physical dynamics.

`PotentialConfig.common` must have no fields set. Classical LAMMPS force
fields do not lower the electronic structure overlay. The selected union arm
must be `lammps`.

## Evaluation

| Function | Contract |
| --- | --- |
| `lammpc_potential_result_size_for_force_input` | Returns safe output capacity without creating LAMMPS. |
| `lammpc_session_calculate_result` | Evaluates on a persistent session. |
| `lammpc_calculate_result` | Creates a temporary session from `LammpsParams`. |
| `lammpc_calculate_result_from_config` | Creates a temporary session from `PotentialConfig`. |
| `lammpc_energy` | Legacy raw array entry point. It fails in current builds. |
| `lammpc_energy_forces` | Legacy raw array entry point. It fails in current builds. |
| `lammpc_energy_gradient` | Legacy raw array entry point. It fails in current builds. |

The session evaluation requires nonnull session, input, result buffer, and
result size output pointers. On a small output buffer, it returns failure and
stores the required capacity in the size output. Geometry requires at least
one atom, exactly three positions per atomic number, and either zero or nine
cell entries.

Cells use row major restricted triclinic form. Entries 1, 2, and 5 must be
zero within `1e-12`. Diagonal entries must be positive. Atomic numbers must
exist in `typeToAtomicNumber`.

## Projected MD

| Function | Contract |
| --- | --- |
| `lammpc_session_md_initialize` | Creates geometry if needed and initializes projected BAOAB state. |
| `lammpc_session_md_step` | Advances exactly `step_count` successful iterations. |
| `lammpc_session_md_state` | Copies selected state arrays and counters. |

Initialization requires nonnull pointers, a coordinate count equal to three
times atom count, a finite nonzero normal, and finite point values. Time step,
thermal energy, and friction must be finite and greater than zero. The seed
has no reserved value. Zero is deterministic.

Mass weighted coordinates interpret `point` in `sqrt(amu) * angstrom` and the
normal in the corresponding reciprocal coordinate direction. Cartesian
coordinates interpret both point and normal in ordinary Cartesian space.

State arrays have `coordinate_count` entries. Each output array and the
observation pointer may be `NULL`. A state request before initialization or
with the wrong count fails.

## NVE

| Function | Contract |
| --- | --- |
| `lammpc_session_nve_set_positions` | Scatters positions and optionally updates a nine entry cell. |
| `lammpc_session_nve_set_velocities` | Scatters Cartesian velocities. |
| `lammpc_session_nve_step` | Runs LAMMPS `fix nve` for the requested count. |
| `lammpc_session_nve_state` | Copies positions, velocities, PE, and KE. |

Every NVE call requires live geometry from a successful evaluation. Array
counts must equal three times the live atom count. A position call without a
cell must pass `NULL` and a cell length of zero. A position call with a cell
must pass exactly nine entries. The NVE time step must be finite and positive,
including when the step count is zero.

## Hard boundaries and direct integration

| Function | Contract |
| --- | --- |
| `lammpc_hard_boundary_apply` | Applies one callback decision to supplied buffers. |
| `lammpc_session_set_hard_boundary` | Stores or clears one session callback. |
| `lammpc_session_clear_hard_boundary` | Clears the callback. |
| `lammpc_session_integrate` | Runs a preconfigured LAMMPS integrator. |

`lammpc_hard_boundary_apply()` returns the callback decision or minus one.
A null callback returns accept without reading any state pointer. Reject
restores both previous arrays. Reflect keeps callback edits. An unknown code
restores both arrays and fails.

Direct integration requires live geometry, a nonnegative step count, and a
LAMMPS time integrator configured through `extraSetup`. Projected MD manages
its own integrator. NVE manages a private `fix nve`.

## Thread safety

The following operations are safe across threads when each thread uses a
different session:

* Session creation, evaluation, dynamics, callbacks, and destruction.
* Reading the immutable feature table.
* Reading each thread's latest error.

No operation is safe concurrently on the same session. A callback runs in the
thread that advances its session. Its `user` pointer must remain valid until
the callback is cleared or the session is destroyed.
