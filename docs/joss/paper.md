---
title: 'ARBFN: Arbitrary Externally-Computed Closed-Loop Force Fixes in LAMMPS'
tags:
  - LAMMPS
  - molecular dynamics
  - MPI
  - HPC
  - simulation
authors:
  - name: Jordan Dehmel
    orcid: 0009-0002-0982-8149
    affiliation: 1
  - name: Jarrod Schiffbauer
    orcid: 0000-0001-9791-421X
    affiliation: "2, 3"
affiliations:
 - name: Dept of Computer Science, Colorado Mesa University, Grand Junction, CO
   index: 1
   ror: 0451s5g67
 - name: Dept of Physical and Environmental Sciences, Colorado Mesa University, Grand Junction, CO
   index: 2
   ror: 0451s5g67
 - name: Dept of Mechanical Engineering, University of Colorado Boulder, Grand Junction, CO
   index: 3
   ror: 02ttsq026
date: 9 May 2025
bibliography: paper.bib
---

# Abstract

The molecular dynamics simulation software LAMMPS
(Large-scale Atomic/Molecular Massively Parallel Simulator)
[@LAMMPS] provides a scripting language for the easy
implementation of experiments: However, it is not
all-encompassing. There are many situations in which LAMMPS
alone is not sufficient and some external computation must
be used (for example, quantum effects [@qmmm]
and machine learning control [@Rohskopf2023]. Previous
researchers, when confronted with these limitations, have
implemented custom interfaces for specific programs. This
project outlines the development of a generic protocol for
this process. Specifically, we introduce an
externally-controlled arbitrary atomic forcing fix within
the existing MPI (Message Passing Interface) LAMMPS
framework [@mpi-docs]. This involves an
arbitrary-language controller program being instantiated
alongside LAMMPS in an MPI runtime, then communicating with
all LAMMPS instances whenever the desired fix must be
computed. The protocol communicates in JSON (JavaScript
Object Notation) strings and assumes no controller
linguistic properties beyond a valid MPI implementation.

# Introduction

LAMMPS fixes are powerful, but sometimes insufficient. Although
modifying LAMMPS' source code allows efficient implementation of
arbitrary fixes, it requires a technical background and imposes
`C++` as the language of choice. Especially when working
with heavyweight or language-specific systems, it would be far
easier to write fixes externally. Thus, we have developed a
system for external force fixes as functions of the entire
simulated system. Furthermore, in cases where such
inter-process-communication-heavy computation would be overkill,
we provide a custom externally-determined forcing field.

If we want to have an external system act as a "controller"
over our LAMMPS particles, we will need to define a
`C++ fix` class which can then be applied in scripts.
Instances of this class will need to be able to communicate
externally: The easiest way to do this is via the Message
Passing Interface (MPI), which LAMMPS already uses. These
"workers" will send MPI packets to the controller whenever
they need an update, then receiving a result and applying it.
For readability and encoding-independence, we will send packets
using JavaScript Object Notation (JSON). To avoid gridlock, we
will allow the user to specify a maximal time to await
controller response before an error is thrown.

We will call this fix type **`fix arbfn`** (for "arbitrary
function" of the state of the simulation).

Note: JSON incurs overhead cost proportional to the size of the
message because of its syntax. It would be faster and smaller to
send raw encodings of the values used, at the cost of imposing
additional restrictions upon the controller language. Thus, we
have chosen to pay the overhead for JSON.

In the aforementioned special cases wherein constant MPI
communication is unnecessary, we will also define the
**`arbfn/ffield`** fix, which determines a static force field
via MPI communication at instantiation, then trilinearly
interpolating atom positions onto a finite position grid in
order to find their forcing values at runtime.

Our final `LAMMPS` interface for `fix arbfn` is
exemplified by the following code.

```
# Every timestep, send all atomic data and
# receive fix data.  If the controller
# takes longer than 50 ms to respond, error.
fix name_1 all arbfn maxdelay 50.0

# Every 100 timesteps, do as above with no
# time limit
fix name_2 all arbfn every 100
```

Likewise, `fix arbfn/ffield` is shown below.

```
# At initialization, retrieve a mesh of 101
# by 201 by 301 nodes. Every timestep,
# perform trilinear interpolation of the
# received force field
fix name_3 all arbfn/ffield 100 200 300

# Every 100 timesteps, send all atom data
# to the controller and refresh the grid
fix name_4 all arbfn/ffield 10 10 10 every 100
```

# Details

## `fix arbfn` \label{arbfn}

The first fix provided by the package is `fix arbfn`. It
is the most powerful and the slowest. Every time this fix is
called, its atoms are sent off to the controller over MPI. The
controller then determines some amount of force to add to each
atom, sending it back to LAMMPS to implement. The controller is
also allowed to send a "waiting" packet indicating that LAMMPS
should wait another few milliseconds. If no response is received
within some specified time limit, LAMMPS will error. Since there
may be arbitrarily many LAMMPS instances running, the controller
may choose to await all the data or to send back data
immediately. It is slightly faster to send back data one
controller at a time, but limits the capabilities of the fix
(for instance, a fix pushing atoms towards the center of mass
could not be implemented).

This fix can be used to implement frame-by-frame control of the
forces of atoms based on the state of some external system. As
a frivolous example, imagine a LAMMPS simulation where forces
could be applied by using a physical joystick. It can also be
used to apply forces based on attributes not feasibly
implementable solely within LAMMPS.

The protocol for `fix arbfn` is shown in figure
\autoref{fig:arbfn_protocol}. Note that communication between LAMMPS
and the "worker" (fix object instance) is virtually free,
while communication between the worker and the controller is
very expensive.

![`fix arbfn` protocol.\label{fig:arbfn_protocol}](arbfn_protocol.png)

While necessary for some use cases, this fix is painfully slow:
A simulation that may take only a few minutes without it will
instead take hours.

## `fix arbfn/ffield` \label{ffield}

Evolving from the aforementioned MPI delays is the
`arbfn/ffield` fix. This takes in some spatial grid of
nodes at instantiation via MPI, then interpolates between them
to find specific force field values.

This fix is *not* able to update frame-by-frame, and the
interpolation it does is position-only (velocity, existing
force, and orientation cannot come into play), but is in
exchange about 100 times faster.

The protocol for the `arbfn/ffield` fix is shown in figure
\autoref{fig:ffield_protocol}. Note that there are no longer
costly MPI calls within the simulation loop, and thus the
simulation will perform much better.

![`fix arbfn/ffield` protocol.\label{fig:ffield_protocol}](ffield_protocol.png)

Although the difference between figures \autoref{fig:arbfn_protocol}
and \autoref{fig:ffield_protocol} may seem trivial, the omission of
the controller from the simulation loop allows
\autoref{fig:ffield_protocol} to run orders of magnitude faster.

## `fix arbfn/ffield` with `every n` \label{ffield_every}

Our final protocol addresses the holes in what
`fix arbfn/ffield` can compute. Instead of receiving a
single interpolation grid at the beginning and using it for the
entire simulation, the `every n` argument allows us to
dynamically update the grid every $n$ time steps. Specifically,
every $n$-th time step, the worker sends all of its atomic data
to the controller and receives a new interpolation grid in
return.

This protocol allows more generality at the cost of speed, while
still being (generally) faster than \autoref{arbfn}. It allows a
degree of dependency of the force field on the atomic data (e.g.
center of mass) which was previously impossible.

The protocol for `fix arbfn/ffield` with the
`every n` argument is shown in figure
\autoref{fig:ffield_every_protocol}. Note that this reintroduces the
costly IPC during the simulation, but is still more sparse than
\autoref{arbfn}.

![`fix arbfn/ffield` protocol when used with the `every n` argument.\label{fig:ffield_every_protocol}](ffield_every_protocol.png)

\autoref{ffield} is a special case of \autoref{ffield_every} where $n$
is larger than the length of the simulation. Internally, this is
represented by `every 0`. The other limit case,
`every 1`, is *nearly* \autoref{arbfn}: It communicates
every frame (and therefore is at least as slow), but the atom
forces are ultimately still interpolated according to the grid,
rather than directly controlled.

## Running Simulations

In order to keep ARBFN MPI communication from interfering with
internal LAMMPS communication, we must run LAMMPS with the
`-mpicolor ...` command-line argument. The number following this
must be anything **except $56789$** (this color is reserved for
internal package communication). On UNIX systems, this takes the
form that follows.

```
mpirun -n 1 ./controller : \
    -n ${num_lmp_threads} \
    lmp -mpicolor 123 \
    -in input_script.lmp
```

# Performance Testing

Our performance experiments were carried out with small
simulations. We kept a constant particle density of 2 atoms per
unit area in a 2D system (where the square box's edge length was
varied) to avoid minimization problems arising from varying
particle count directly. The controllers used were minimal so as
to test LAMMPS performance rather than external computation
time: Although they always applied at least some forces to some
atoms, no significant computation was involved. These
experiments are intended to demonstrate the trend of the running
times in small simulations, as well as to corroborate that our
fixes scale proportionally to a no-fix system. A `python`
script was used to automate the running of scripts.

![Comparing `fix arbfn` to control and `fix arbfn/ffield`. While still appearing linear with respect to the number of atoms, it runs *much* slower than the latter two (which are about the same speed).\label{fig:arbfn_comparison}](arbfn_scaled_comparison.png)

Figure \autoref{fig:arbfn_comparison} demonstrates the sharp slope
of `fix arbfn`: IPC is extremely costly, and
frame-by-frame updates should be avoided whenever possible. That
being said, the time usage appears to scale proportionally to
control as expected. The `fix arbfn/ffield` results
appear to scale with a much smaller coefficient, as expected:
This is a much more usable fix, although its use case is more
narrow.

# Conclusion

We discussed the implementation of the `ARBFN` package for
LAMMPS, including a brief analysis of its performance as
simulations scaled. This package allows LAMMPS fixes which are
determined at runtime by arbitrary external "controllers",
either one-and-done (`arbfn/ffield`) or frame-by-frame
(`fix arbfn`). Initial testing shows that the former
performs nearly as well as unmodified LAMMPS, while the
latter performs about 2 orders of magnitude worse. We also
provided an option to update the interpolation force field
as a function of atom data periodically throughout the
simulation's life cycle. The package is provided as
supplementary material with this paper.

# Related Work

The source code of this package draws from the QMMM package
[@qmmm] and the LAMMPS codebase [@LAMMPS] for the
framework of MPI communication. It also used the
`BROWNIAN` package, "Extending and Modifying LAMMPS"
[@modifying-lammps], and [@chemengineering5020030] as a
basis for force modifications. Packages like QMMM
and FitSNAP [@Rohskopf2023] use custom `C++` to
interface with external controllers: Our package seeks to
standardize such communication. The external interfacing of our
package is similar to the built-in LAMMPS `python`
wrapper (used to similar effect in e.g. [@do2024feedback]),
but allows more linguistic generality. Indeed, `python`
has proved a popular choice in the composure and control of MD
simulations: Source [@balasubramanian2016extasy] uses
template `python` scripts to control and scale
simulations.

Visualization has been hand-in-hand with control as a primary
concern in molecular dynamics. Systems like
[@koutek2002virtual] and [@stone2001system] allow
hand-modification of particle forces within their systems, even
though their concern is primarily in human-simulation interface
rather than software-simulation.

# Acknowledgements

Supported by NSF grant 2126451. Based upon work partially
funded by Colorado Mesa University.

# References
