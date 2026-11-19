# Design explanation

## Why sessions exist

LAMMPS owns mutable simulation state. Creating that state for every force call
would repeat box, atom, mass, and pair style setup. A `LammpcSession` keeps one
LAMMPS instance and one decoded interaction model. Geometry updates can then
write positions into the existing instance.

The C type is opaque so its C++ layout does not become part of the ABI. The
host owns the pointer and ends its lifetime with
`lammpc_session_destroy()`. Reconfiguration destroys the live LAMMPS instance
because pair state and atom types can depend on every configuration field.

## Wire and unit boundaries

The wire format comes from the pinned `Potentials.capnp`. Configuration and
evaluation messages are flat Cap'n Proto buffers. lammpc converts input
geometry to angstroms, evaluates a LAMMPS instance in `metal` units, and
converts the result to the units requested by `ForceInput`.

LAMMPS unit styles change more than labels. They define constants, time units,
and accepted pair data. The implementation therefore rejects `real` instead
of applying partial conversions. `metal` is the only supported style.

## Projected BAOAB

The constrained stepper uses the BAOAB splitting for Langevin dynamics.
Leimkuhler and Matthews analyze this splitting and its configurational
sampling properties in
[DOI 10.1063/1.4802990](https://doi.org/10.1063/1.4802990).

For a mass weighted plane, lammpc defines

```text
q_i = sqrt(m_i) x_i
n dot (q - point) = 0
```

The supplied normal is normalized in Euclidean `q` coordinates. Position,
velocity, force, and random increments are projected. For a Cartesian plane,
the constraint is `n dot (x - point) = 0` and projection uses the inverse mass
metric.

Initialization draws Maxwell velocities from the configured seed and thermal
energy. Equal seeds and equal inputs produce byte identical trajectories on
the same build. A rejected hard boundary proposal restores the physical state
but retains consumed random draws. This makes the next proposal use fresh
noise. An invalid callback return restores the random state because the step
failed.

## NVE propagation

The NVE API delegates propagation to LAMMPS `fix nve`, which uses velocity
Verlet. Thompson and coauthors describe the current LAMMPS architecture in
[DOI 10.1016/j.cpc.2021.108171](https://doi.org/10.1016/j.cpc.2021.108171).

lammpc installs one private NVE fix on first use and updates the LAMMPS time
step for every call. State reads synchronize positions, velocities, potential
energy, and kinetic energy from the live instance.

Periodic coordinates may arrive several cells away. lammpc maps periodic
fractional components closer to the live cell before scattering atoms.
Nonperiodic components remain unchanged. The mapping uses the full restricted
triclinic lattice, including tilt factors.

## Hard boundary ownership

lammpc transports state to a host callback but does not define region
geometry. Accept keeps the post step state. Reject restores the previous
state. Reflect keeps callback edits and, for projected dynamics, restores the
plane constraint before evaluating a new force.

The callback belongs to one session. The direct LAMMPS integration path
attaches a temporary LAMMPS fix only while `lammpc_session_integrate()` runs.
The projected stepper invokes the same callback from its own loop. NVE does not
invoke it.

## Concurrency model

The diagnostic buffer is thread local in the embed (a 512 byte buffer updated
only by session creation, configuration, evaluation, MD, NVE, hard boundary
hook, and direct integration calls) and a read only static string in the stub.
Hard boundary state and all dynamics state belong to a session. This permits
independent sessions in independent threads.

A session has no internal mutex. Concurrent access to one session can race on
its LAMMPS handle, vectors, counters, and callback state. External
synchronization is required when a session pointer can reach more than one
thread. Process wide LAMMPS limitations may also depend on the packages linked
into the chosen LAMMPS build.
