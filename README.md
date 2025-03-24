
# Arbitrary Externally-Computed Force Fixes in LAMMPS (`ARBFN`)

![Test Badge](https://github.com/jorbDehmel/lammps_arb_fn/actions/workflows/ci-test.yml/badge.svg)

J Dehmel, J Schiffbauer, Colorado Mesa University 2024/25 (MIT License)

## Motivation and Abstract

When dielectric colloidal Janus particles are placed in
pseudo-2D environments, severe wall effects come into play
complicating the measurement of bulk motility. The efficient and
easy implementation of semi-empirical wall-effect systems in
LAMMPS would ease the understanding of these effects and
possibly lead to their mitigation or control. In this case, the
addition of arbitrary force fields onto the simulation space
would be useful. This project introduces an
externally-controlled arbitrary atomic forcing fix for LAMMPS
within the existing MPI framework. This involves an
arbitrary-language controller program being instantiated
alongside LAMMPS in an MPI runtime, then communicating with all
LAMMPS instances whenever the desired fix must be computed. The
protocol communicates in JSON strings and assumes no linguistic
properties except a valid MPI implementation.

## Requirements and Testing

This software is built to work with the LAMMPS source code, and
is such written in `C++` using MPI. The following are
requirements for compilation:

- `g++` and `make` (you almost certainly already have these)
- `CMake`
- Some MPI provider (EG `openmpi`: You probably also have this
    if you are compiling LAMMPS)
- `boost` library for `C++` (specifically `boost_json`)
- `python3` (for the install script)

The following software is required for testing, but not
necessarily for non-testing compilation.

- `python3-pip`
- `mpi4py`
- `git`
- (Optional) `docker` or `podman`

To check your system for the given requirements,
run: `make check`

To test the source code, run: `make test`

(More details on installing the package inside a LAMMPS binary
will be given later)

To enter a development `Docker` container, run `make docker`

To enter a development `podman` container, run `make podman`

After compiling and install LAMMPS with the extension, a simple
testing script can be found in `./lmp_test`. You can run it
(on a 4-core or more system) via `make run`.

## Installation

This package is currently only provided standalone, and thus
requires **a full re-compilation of LAMMPS from source**. This
document will assume a `CMake` build, rather than a traditional
`make` build. Please see
[the LAMMPS docs on building from source](https://docs.lammps.org/Build.html)
for more details.

Since this process is not already in the LAMMPS `CMake` system,
it provides `INSTALL.py` locally to facilitate installation.
This script first locates the LAMMPS source code repository,
then copies the `ARBFN` source code into it, then modifies
`CMake` so that it will build properly, then optionally runs the
`CMake` installation process.

To run the installation script, simply run

```sh
python3 INSTALL.py
```

from this directory and follow the instructions given.

**Note:** this requires you to have already cloned the LAMMPS
source code somewhere locally: Directions for this can be found
in the aforementioned LAMMPS build guide.

## Running Simulations

Although LAMMPS is built on MPI, extra care is needed when
running it alongside other processes. Namely, the `-mpicolor X`
command-line argument must be used, where `X` is any number
**except** 56789 (with 56789 being reserved for internal ARBFN
communication). Note that
**`-mpicolor` must be the first argument after `lmp`**.
Additionally, `lmp` must be launched by an `mpirun` call
**wherein the controller is also launched**. Only one instance
of the controller should be used, but arbitrarily many instances
of LAMMPS are allowed.

For instance, to launch one instance of the local controller
executable `./controller` alongside 3 instances of LAMMPS with
color $123$ on `input_script.lmp`:

```sh
mpirun -n 1 \
    ./controller \
    : -n 3 \
    lmp -mpicolor 123 -in input_script.lmp
```

To use the same setup, but instead use the `python` controller
`./py_controller.py`:

```sh
mpirun -n 1 \
    python3 ./py_controller.py \
    : -n 3 \
    lmp -mpicolor 123 -in input_script.lmp
```

## Using the `arbfn` Fix

The `ARBFN` package provides the `arbfn` fix, shown below.

```lammps
fix name_1 all arbfn
fix name_2 all arbfn maxdelay 50.0
fix name_3 all arbfn every 100
fix name_4 all arbfn maxdelay 50.0 every 100
```

The `maxdelay X` (where `X` is the max number of milliseconds to
await a response before erroring, with $0.0$ ms meaning no
limit) and `every Y` (where the fix is applied every `Y`
timesteps, with $1$ being every step and $0$ being undefined)
arguments are both optional. The default max delay is $0.0$ (no
limit) and the default periodicity is $1$ (apply every
time step).

There is also the `dipole` argument, which includes the values
`"mu"`, `"mux"`, `"muy"`, `"muz"` from LAMMPS for each atom.

```lammps
fix name_5 all arbfn dipole
```

## `fix arbfn` Protocol

This section uses pseudocode and standard MPI calls to outline
the protocol used from both the controller and worker
perspective. The reader should assume that there is exactly one
controller and at least one LAMMPS "worker".

We begin by describing the protocol from the controller's
perspective.

1. Initial MPI setup
    - Call `MPI_Init` to initialize the system
    - Call `MPI_Comm_split` to split `MPI_COMM_WORLD` off into
        a "junk" communicator which can then be discarded. This
        is necessary because LAMMPS uses an internal
        `MPI_Comm_split` upon instantiation, and all processes
        in the world must make the call before the process will
        advance. If this is not done, the system will hang
        without error indefinitely.
    - Call `MPI_Comm_split` **for the 2nd time**, this time
        splitting `MPI_COMM_WORLD` off using the color $56789$
        (the color all our MPI comms are expected to have) into
        a communicator which we save. This will be the
        communicator that we use for the remainder of the
        protocol, and corresponds to the splitting off of the
        `ARBFN` fixes from the default LAMMPS communicator. The
        same synchronization issues will occur upon omission of
        this step as the previous.
2. Enter server loop
    - Unless exited, repeat this step (2) forever after completion
    - Make a call to `MPI_Probe` with any source and tag,
        storing the resulting information. Since we are using
        the colored communicator, this will only every receive
        messages from the worker instances.
    - Upon receiving the non-null-terminated ASCII string
        message, save it to a string and decode it into a JSON
        object. The JSON object is guaranteed to have an
        attribute with the key `"type"`.
        - If `"type"` is the string `"register"`, increment some
            counter of the number of registered workers and send
            back a JSON packet with type `"ack"`.
        - If `"type"` is the string `"deregister"`, decrement
            the aforementioned counter. If it is now zero, exit
            the server loop. This is the only case in which the
            server shuts down.
        - If `"type"` is the string `"request"`, the JSON will
            encode the data (at minimum "x", "y", "z", "vx",
            "vy", "vz", "fx", "fy", and "fz": In some cases also
            the dipole information "mux", "muy", "muz") of each
            atom it owns into a list with the key `"atoms"`. The
            controller is expected to respond with either a
            packet of type `"waiting"` (requiring no additional
            information but prompting the worker to resent the
            packet after some interval) or a packet of type
            `"response"`. A `"response"` packet will have a list
            called `"atoms"` where the $i^\texttt{th}$ entry
            corresponds to the $i^\texttt{th}$ atom in the
            request. Each atom in the response will have (at
            minimum) "dfx", "dfy", and "dfz". Each of these will
            be a double corresponding to the prescribed deltas
            in force for their respective dimension.
3. Shutdown
    - After all workers have send `"deregister"` packets, LAMMPS
        will begin shutting down. This entails one final MPI
        barrier, so we the controller must call `MPI_Barrier`
        on **`MPI_COMM_WORLD`**. After this, LAMMPS will halt.
    - Now, the controller must shut down and free its resources
        in the traditional `C` MPI way (EG `MPI_Comm_free` and
        `MPI_Finalize`). The server can now halt. If improperly
        synced, the system will hang without error.

This ends our description of the controller protocol. We will
now describe the worker side of the protocol, omitting any
information which can be deduced from the above.

1. LAMMPS internal setup
    - Since our workers are fix objects, they have no say in the
        initial MPI setup of LAMMPS. This is where the first
        `MPI_Comm_split` call synchronization occurs.
2. Additional setup
    - Upon fix object instantiation, we must make a second
        `MPI_Comm_split` splitting `MPI_COMM_WORLD` into a
        usable communicator with the color $56789$. This
        corresponds with our second synchronization call on the
        controller side.
3. Controller discovery
    - At instantiation, our fix does not know the rank of the
        controller. Thus, the worker will iterate through all
        non-self ranks in the communicator and send a
        `"registration"` packet to each one. All but one of the
        recipients will be workers and not respond, but the one
        that sends back an `"ack"` packet will be recorded as
        the controller.
4. Work
    - For as long as LAMMPS lives, it will call the fix to do
        work in the form of the `post_force` procedure. This
        will iterate through the owned atoms of this worker
        instance, send them in the aforementioned format to the
        controller, receive the controller's prescription, and
        add the force deltas.
5. Deregistration and cleanup
    - Once LAMMPS is done, the fix destructor will be called.
        This method must send the deregistration packet to the
        controller and free up any resources used (MPI or
        standard). The worker **does not** need to call
        `MPI_Barrier`, unlike the controller.

When developing a controller, it is best to use the provided
example controllers in `./tests/` as templates.

## Using the `arbfn/ffield` Fix

The `ARBFN` package also provides the `arbfn/ffield` fix, shown
below.

```lammps
fix name_1 all arbfn/ffield 100 200 400
fix name_3 all arbfn/ffield 1 2 3 every 100
```

Where the number of $x$ **bins** is $100$ (meaning $101$ nodes),
the number of $y$ bins is $200$, and the number of $z$ bins is
$400$. `every n` means that the worker will pause every $n$
time steps to send all atomic data to the controller and refresh
its interpolation grid. The default `every` value is $0$,
meaning that a new grid is never requested: The first grid is
always used instead.

## `fix arbfn/ffield` Protocol

This section uses pseudocode and standard MPI calls to outline
the protocol used from both the controller and worker
perspective. The reader should assume that there is exactly one
controller and at least one LAMMPS "worker".

We begin by describing the protocol from the controller's
perspective.

1. Initial MPI setup
    - Call `MPI_Init` to initialize the system
    - Call `MPI_Comm_split` to split `MPI_COMM_WORLD` off into
        a "junk" communicator which can then be discarded. This
        is necessary because LAMMPS uses an internal
        `MPI_Comm_split` upon instantiation, and all processes
        in the world must make the call before the process will
        advance. If this is not done, the system will hang
        without error indefinitely.
    - Call `MPI_Comm_split` **for the 2nd time**, this time
        splitting `MPI_COMM_WORLD` off using the color $56789$
        (the color all our MPI comms are expected to have) into
        a communicator which we save. This will be the
        communicator that we use for the remainder of the
        protocol, and corresponds to the splitting off of the
        `ARBFN` fixes from the default LAMMPS communicator. The
        same synchronization issues will occur upon omission of
        this step as the previous.
2. Enter server loop
    - Unless exited, repeat this step (2) forever after
        completion
    - Make a call to `MPI_Probe` with any source and tag,
        storing the resulting information. Since we are using
        the colored communicator, this will only every receive
        messages from the worker instances.
    - Upon receiving the non-null-terminated ASCII string
        message, save it to a string and decode it into a JSON
        object. The JSON object is guaranteed to have an
        attribute with the key `"type"`.
        - If `"type"` is the string `"register"`, increment some
            counter of the number of registered workers and send
            back a JSON packet with type `"ack"`.
        - If `"type"` is the string `"deregister"`, decrement
            the aforementioned counter. If it is now zero, exit
            the server loop. This is the only case in which the
            server shuts down.
        - If `"type"` is the string `"gridRequest"`, the JSON
            will have the values `"offset"` (the x/y/z box
            starting coordinates), `"spacing"` (the dx/dy/dz
            between nodes in space), and `"nodeCounts"` (the
            number of nodes per side). It will expect a return
            value of an object with `"nodes"` being an array.
            Each item of `"points"` must have attributes
            `xIndex`, `yIndex`, `zIndex` (the bin indices),
            `dfx`, `dfy`, and `dfz` (the force field
            contributions at that node). If this is not the
            initial grid request, there will also be an
            `"atoms"` attribute to the request, just as in
            `fix arbfn`.
3. Shutdown
    - After all workers have send `"deregister"` packets, LAMMPS
        will begin shutting down. This entails one final MPI
        barrier, so we the controller must call `MPI_Barrier`
        on **`MPI_COMM_WORLD`**. After this, LAMMPS will halt.
    - Now, the controller must shut down and free its resources
        in the traditional `C` MPI way (EG `MPI_Comm_free` and
        `MPI_Finalize`). The server can now halt. If improperly
        synced, the system will hang without error.

This ends our description of the controller protocol. We will
now describe the worker side of the protocol, omitting any
information which can be deduced from the above.

1. LAMMPS internal setup
    - Since our workers are fix objects, they have no say in the
        initial MPI setup of LAMMPS. This is where the first
        `MPI_Comm_split` call synchronization occurs.
2. Additional setup
    - Upon fix object instantiation, we must make a second
        `MPI_Comm_split` splitting `MPI_COMM_WORLD` into a
        usable communicator with the color $56789$. This
        corresponds with our second synchronization call on the
        controller side.
3. Controller discovery
    - At instantiation, our fix does not know the rank of the
        controller. Thus, the worker will iterate through all
        non-self ranks in the communicator and send a
        `"registration"` packet to each one. All but one of the
        recipients will be workers and not respond, but the one
        that sends back an `"ack"` packet will be recorded as
        the controller.
4. Initialization and work
    - At instantiation, the worker will send a `"gridRequest"`
        packet to the server. This will contain the data
        mentioned in the previous section, and the controller
        will respond in the aforementioned way
    - After getting the data for each node, the worker will save
        it locally
    - When an atom needs fixed, the worker will find the 8
        nearest mesh points. It will perform a trilinear
        interpolation in $(x, y, z)$ space for `dfx`, `dfy`, and
        `dfz`, adding the results to the atom's total forces. If
        `every n` was specified at instantiation, we check if
        the current time step is a multiple of $n$ **first**. If
        it is, we refresh the interpolation grid
        **before updating.**
5. Deregistration and cleanup
    - Once LAMMPS is done, the fix destructor will be called.
        This method must send the deregistration packet to the
        controller and free up any resources used (MPI or
        standard). The worker **does not** need to call
        `MPI_Barrier`, unlike the controller.

When developing a controller, it is best to use the provided
example controllers in `./tests/` (for `Python` and `C++`) as
templates. These contain built-in timeouts to handle cases where
LAMMPS enters an undetected error state, and therefore avoids
eating away valuable server time.

## Disclaimer

FOSS under the MIT license. Supported by NSF grant 2126451.
Based upon work partially funded by Colorado Mesa University.
