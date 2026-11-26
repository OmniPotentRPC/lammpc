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
mole, and kilocalories per mole. An empty value means electronvolt (`eV`).
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
failure. `lammpc_hard_boundary_apply()` returns a decision code on success
(accept, reject, or reflect) or minus one on error. `lammpc_abi_version()` and
`lammpc_available()` return numeric values directly.

In the embed build, `lammpc_last_error()` returns a nonnull, null terminated
pointer into a 512 byte thread local buffer (messages longer than 511 bytes are
truncated). Its contents are overwritten or cleared by the next error updating
call on the same thread, and the buffer lifetime ends when that thread exits.
Copy the string if needed before another updating call or thread exit, and do
not free or modify it.

The functions that set `lammpc_last_error()` on failure and clear it to `""` on
success in the embed build are `lammpc_configure()`,
`lammpc_session_create()`, `lammpc_session_create_from_config()`,
`lammpc_session_set_params()`, `lammpc_session_configure()`,
`lammpc_session_calculate_result()`, `lammpc_calculate_result()`,
`lammpc_calculate_result_from_config()`, all `lammpc_session_md_*()` and
`lammpc_session_nve_*()` functions, `lammpc_session_set_hard_boundary()`,
`lammpc_session_clear_hard_boundary()`, and `lammpc_session_integrate()`.
`lammpc_set_params()` and the legacy raw array functions (`lammpc_energy()`,
`lammpc_energy_forces()`, `lammpc_energy_gradient()`) always fail. In the embed
build they store an error. In the stub build `lammpc_set_params()` returns `-1`
and the raw array functions return `ok = 0` with the static unavailable
message in `LammpcResult.message`; neither sets `lammpc_last_error()`.

The following calls neither set nor clear `lammpc_last_error()`:

| Function | Reason |
| --- | --- |
| `lammpc_hard_boundary_apply` | Pure C callback helper; returns decision code or `-1`. |
| `lammpc_capabilities_result` | Pure C encoder. Returns `-1` and writes nothing when `capabilities_capnp_size_bytes` is `NULL`. Otherwise it sets `*capabilities_capnp_size_bytes` to zero first. If allocating or encoding the message fails, including a `capn_set_text` or `capn_setp` failure, it returns `-1` and the size stays zero. If encoding succeeds, it stores the required byte count there, then returns `-1` when the output buffer is `NULL` or smaller than that count, and otherwise copies the message and returns `0`. |
| `lammpc_abi_version` | Returns `LAMMPC_ABI_VERSION` directly. |
| `lammpc_available` | Returns `1` or `0` directly. |
| `lammpc_potential_result_size_for_force_input` | Returns byte count or `0` on invalid or unrepresentable input in the embed, and always returns `0` in the stub. |
| `lammpc_session_destroy` | Releases session state; returns `void`. |
| `lammpc_finalize` | No operation; returns `void`. |
| `lammpc_version` | Returns static version text. |
| `lammpc_last_error` | Reads the current diagnostic buffer without modifying it. |
| `lammpc_feature_count`, `lammpc_feature_table`, `lammpc_feature_find` | Read-only static feature table queries. |

In the stub build, `lammpc_last_error()` always returns a fixed read only
static unavailable string, and no call sets or clears it.

Evaluation functions return `LammpcResult`. `ok` is one on success and zero
on failure. `message` is always null terminated. On success, it contains
`ok` and the embed clears `lammpc_last_error()`. On failure, `message` and
`lammpc_last_error()` both carry the diagnostic.

Creation functions return `NULL` on failure. The result size query returns
zero for an invalid `ForceInput`, an unrepresentable size, or any call on the
stub build.

No C function throws a C++ exception across the ABI. As with ordinary C,
passing a dangling pointer or a buffer smaller than its declared size is
outside the contract.

## Lifecycle and discovery

| Function | Contract |
| --- | --- |
| `lammpc_version` | Returns library owned static text. The embed string is `lammpc/` plus the project version. The stub string is `lammpc-stub/` plus that version. |
| `lammpc_abi_version` | Returns `LAMMPC_ABI_VERSION`, currently zero. |
| `lammpc_available` | Returns one for the embed and zero for the stub. |
| `lammpc_last_error` | Returns a nonnull pointer to the calling thread's 512 byte error buffer in the embed, or a fixed static string in the stub. |
| `lammpc_finalize` | No operation in the current ABI. |
| `lammpc_capabilities_result` | Writes the backend `Capabilities` message. |
| `lammpc_feature_count` | Returns the number of static feature entries. |
| `lammpc_feature_table` | Returns a library owned static table. |
| `lammpc_feature_find` | Returns a matching static entry or `NULL`. |

For a capabilities size query, pass `NULL`, zero, and a valid size output
pointer. The function returns minus one and writes the required byte count.
Repeat the call with a buffer of at least that size. A missing size output
pointer is an error. If allocating or encoding the message fails, the function
also returns minus one, but the size output stays zero.

Feature table pointers and their strings remain valid for the life of the
loaded library. Callers must not free or modify them.

## Configuration and session lifetime

| Function | Contract |
| --- | --- |
| `lammpc_set_params` | Legacy global entry point. Both builds ignore both arguments and return `-1`. The embed stores a message saying that global parameters are not supported in `lammpc_last_error()`. The stub does not modify it. |
| `lammpc_configure` | Validates a `PotentialConfig` by creating and destroying a temporary session. |
| `lammpc_session_create` | Creates an idle session from `LammpsParams`. |
| `lammpc_session_create_from_config` | Creates an idle session from the `lammps` arm. |
| `lammpc_session_set_params` | Replaces parameters and destroys live geometry. |
| `lammpc_session_configure` | Replaces configuration and destroys live geometry. |
| `lammpc_session_destroy` | Releases callback, LAMMPS, and session state. `NULL` is accepted. |

The parameter defaults are `metal`, `atomic`, and `p p p`. `pairStyle` and
`typeToAtomicNumber` are required. An omitted `masses` entry uses the
standard atomic weight of that type's atomic number. Atomic numbers outside
1 through 118 are rejected. An explicit mass replaces the standard weight.

`PotentialConfig.common` must have no fields set. Classical LAMMPS force
fields do not lower the electronic structure overlay. The selected union arm
must be `lammps`.

## Evaluation

| Function | Contract |
| --- | --- |
| `lammpc_potential_result_size_for_force_input` | Returns safe output capacity without creating LAMMPS in the embed, or `0` in the stub. Does not set or clear `lammpc_last_error()`. |
| `lammpc_session_calculate_result` | Evaluates on a persistent session. Writes `PotentialResult` size output as described below. |
| `lammpc_calculate_result` | Creates a temporary session from `LammpsParams`, then evaluates with the same size output rules. |
| `lammpc_calculate_result_from_config` | Creates a temporary session from `PotentialConfig`, then evaluates with the same size output rules. |
| `lammpc_energy` | Legacy raw array entry point. It fails in current builds. |
| `lammpc_energy_forces` | Legacy raw array entry point. It fails in current builds. |
| `lammpc_energy_gradient` | Legacy raw array entry point. It fails in current builds. |

The evaluation entry points require a nonnull session (except the one-shot
helpers, which create a temporary session first), valid `ForceInput` bytes, and
a nonnull `potential_result_capnp_size_bytes` pointer. The `PotentialResult`
buffer may be `NULL` when probing capacity once the `ForceInput` has been
validated as described below.

The size output follows one sequence on all three calls. A null session (on the
session call) or a null size pointer returns failure with `invalid arguments`
and leaves `*potential_result_capnp_size_bytes` unchanged. The one-shot calls
check the size pointer before they create the temporary session, so a null size
pointer never starts session creation. Once the arguments pass, the size output
is set to zero, and for the one-shot calls that happens before session
creation. If the session cannot be created, the call returns that error with
the size output at zero.

`ForceInput` has no separate header to parse. The required flat byte count is
written to the size output only after the message has been read, its atom
count, position count and cell entries have passed the geometry checks, and
both `lengthUnit` and `energyUnit` have resolved to supported units. A message
that cannot be read, fails the geometry checks, or names an unsupported
`lengthUnit` or `energyUnit` therefore returns failure with the size output
still at zero.

After the required count is stored, a `NULL` result buffer or a
`potential_result_capnp_capacity_bytes` smaller than that count returns failure
with `PotentialResult buffer too small` and keeps the required count in the
size output. A failure while copying the geometry, building or evaluating the
LAMMPS instance, or encoding the result also keeps the required count. The
encoding step returns `PotentialResult write failed` when allocating or
encoding the message fails, and it leaves the size output at the required count
because the count written replaces it only on success. On success the call
copies the encoded message and sets the size output to the byte count written.

`lammpc_potential_result_size_for_force_input()` is the usual way to size a
buffer before allocation because it parses geometry and units without a session or
LAMMPS instance. The size output on the evaluation calls still reports an
undersized buffer at call time. In the stub build, evaluation calls fail with
the fixed unavailable message and do not modify the size output.

Geometry requires at least one atom, exactly three positions per atomic
number, and either zero or nine cell entries.

Cells use row major restricted triclinic form. Every entry must be finite.
Entries 1, 2, and 5 must be zero within `1e-12`. Diagonal entries must be
positive. Atomic numbers must
exist in `typeToAtomicNumber`. A number listed once is shared by every atom
of that number. A number listed on several types is assigned in list order,
one unused type per atom, and another atom after those types are used is
rejected.

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
