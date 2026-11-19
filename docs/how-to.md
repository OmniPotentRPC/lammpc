# How to use sessions

## Configure a session

Encode `LammpsParams` or a `PotentialConfig` with the pinned
`Potentials.capnp` schema. A minimal `LammpsParams` message supplies:

1. `pairStyle`.
2. At least one entry in `typeToAtomicNumber`.
3. Any `pairCoeffs` required by the selected pair style.
4. Positive `masses` when molecular dynamics will use values other than the
   one atomic mass unit fallback.

Leave `unitsStyle` empty or set it to `metal`. Leave `atomStyle` empty or set
it to `atomic`. The embed rejects every other value, including
`unitsStyle = "real"`.

Create one session for each independent LAMMPS instance:

```c
LammpcSession *session =
    lammpc_session_create_from_config(config_bytes, config_size);
if (session == NULL) {
  fprintf(stderr, "%s\n", lammpc_last_error());
  return 1;
}
```

Always release it with `lammpc_session_destroy(session)`. Destroying `NULL` is
safe.

## Evaluate energy and forces

Encode a `ForceInput` with one atomic number per atom and three coordinates per
atom. The optional cell has nine row major entries. Ask for an output capacity
before allocating the result:

```c
size_t capacity =
    lammpc_potential_result_size_for_force_input(input, input_size);
void *output = malloc(capacity);
if (output == NULL)
  return 1;
size_t output_size = 0;
LammpcResult result = lammpc_session_calculate_result(
    session, input, input_size, output, capacity, &output_size);
if (!result.ok) {
  fprintf(stderr, "%s\n", result.message);
  free(output);
  return 1;
}
```

The first evaluation creates atoms and the LAMMPS box. Later evaluations reuse
the live instance when atom count and cell match. `PotentialResult` uses the
energy and length units requested by `ForceInput`. The scalar
`LammpcResult.energy_h` always uses Hartree.

The complete executable version of this flow is
`tests/test_embed_lj_e2e.c`.

## Run projected Langevin dynamics

First create a configured session. Then call
`lammpc_session_md_initialize()` with a valid `ForceInput`, a plane point, a
nonzero plane normal, and one coordinate for each Cartesian degree of
freedom. Every control in `LammpcMdConfig` must be finite and positive except
the seed, which may be zero.

Use picoseconds for time, electronvolts for thermal energy, reciprocal
picoseconds for friction, angstroms for Cartesian positions, and atomic mass
units for masses. The stepper is BAOAB with force, velocity, and stochastic
increments projected onto the selected plane.

```c
LammpcMdConfig md = {
    .time_step_ps = 0.0002,
    .thermal_energy_ev = 0.01,
    .friction_per_ps = 10.0,
    .seed = 17,
    .coordinates = LAMMPC_PLANE_MASS_WEIGHTED,
};
if (lammpc_session_md_initialize(
        session, input, input_size, point, normal, coordinate_count, &md) != 0) {
  fprintf(stderr, "%s\n", lammpc_last_error());
  return 1;
}
```

Call `lammpc_session_md_step()` to advance an exact number of iterations.
Call `lammpc_session_md_state()` to copy any required arrays. Output pointers
may be `NULL`, but `coordinate_count` must still match the initialized state.

## Run unconstrained NVE dynamics

An NVE session needs live geometry. Establish it with
`lammpc_session_calculate_result()` before calling an NVE function.

1. Call `lammpc_session_nve_set_positions()`.
2. Call `lammpc_session_nve_set_velocities()`.
3. Call `lammpc_session_nve_step()`.
4. Read the state with `lammpc_session_nve_state()`.

All arrays contain `3 * atom_count` doubles. Positions use angstroms,
velocities use angstroms per picosecond, the time step uses picoseconds, and
observed energies use electronvolts. A zero step count performs a LAMMPS
`run 0`. The time step must still be finite and positive.

See `tests/test_embed_nve_e2e.c` for an EAM trajectory, an energy check, and a
time reversal check.

## Add a hard boundary

Call `lammpc_session_set_hard_boundary()` with a function pointer and an
optional user pointer. The callback returns accept, reject, or reflect. A
reflection means that the callback has already edited both state buffers.
Any other return code is an error.

The callback is active for `lammpc_session_integrate()` and projected MD.
NVE entry points ignore it. Clear it before releasing callback state:

```c
if (lammpc_session_clear_hard_boundary(session) != 0) {
  fprintf(stderr, "%s\n", lammpc_last_error());
  return 1;
}
```

See [Hard boundary ABI](hard-boundary.org) for the decision rules.

## Use sessions from threads

Use a different session in each thread. Independent sessions, callbacks, and
thread local error messages can run concurrently. Never call two functions on
the same session at the same time. The library does not lock session state.
Do not destroy a session while another thread is using it.

The ThreadSanitizer lane exercises four independent sessions. Run it with:

```sh
pixi run -e embed test-embed-tsan
```
