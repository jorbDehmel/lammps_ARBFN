
# The `ARBFN` Package for LAMMPS
## Runtime Arbitrary External Forcing Functions

J Dehmel, Colorado Mesa University, 2024/25

---

# Motivation

Though expressive, LAMMPS forcing functions have limits.

- Might require external computation
- Might want to periodically dump live data to a general-purpose
    language for processing
- Might just be easier to write in another language

**Solution:** Allow outsourcing of force computations to an
external "controller"

---

# Restrictions and Solutions

- Controller should be as arbitrary as possible
    - **Solution**: Define only protocol, not controller
- Protocol should support cross-platform communication
    - **Solution**: JSON packets over MPI
- Must handle latency to avoid wasting server time
    - **Solution**: Throw error if controller isn't timely
- Must assume that LAMMPS is distributed across many nodes
    - **Solution**: Use MPI, which LAMMPS already has
- Must not interfere with LAMMPS' internal communication
    - **Solution**: Use MPI colors to separate communication
        channels

---

# Implementing a "Worker" Fix

- When the fix is instantiated:
    - Send a registration packet to the controller

- When LAMMPS calls our fix:
    - Collect atom data
    - Send a request w/ the atom data to the controller
    - Await a response, erroring if it takes too long
    - Parse the response to get the force deltas
    - Implement the force changes received

- When LAMMPS shuts down:
    - Send a deregistration packet to the controller

---

# Implementing a "Controller" Program

- When the controller is starting:
    - Until all workers have deregistered, listen for packets

- When a worker registers:
    - Log them and send an acknowledgement packet in response

- When a worker deregisters:
    - Log that they have left
    - If this is the last one, halt

- When a worker requests a fix:
    - Receive its atom data
    - If we need the rest of the workers' data, wait for that
    - If we have or don't need the other data, calculate and
        send the force deltas

---

# Using `fix arbfn`

1. Compile LAMMPS with this package enabled (more details in the
    repo root dir)
2. In your LAMMPS script, apply the fix:

```lammps
# Let all atoms be adjusted by the controller every 5 timesteps,
# with a max delay of 50ms.
fix fix_name all arbfn maxdelay 50.0 every 5
```

3. Run the script and controller in the same `mpirun` command

```sh
# "123" can be any number EXCEPT 56789
mpirun -n 1 ../example_controller.out : \
    -n 3 lmp -mpicolor 123 -in test.lmp
```

---

# Example pt. 1: Controller (`wall_effect.cpp`)

```cpp
#include "controller.hpp"
#include <boost/json/object.hpp>
#include <cmath>
constexpr double sigmoid(const double _x) {
  return 1.0 / (1.0 + exp(-_x));
}
int main() {
  // For each atom (as JSON), get the dfx, dfy, and dfz
  independent_controller([](const auto &_json,
                            double &, double &_dfy,
                            double &) {
    auto y = _json.at("y").as_double();
    _dfy = 0.1 * (sigmoid((5.0 - 30.0 - y)) -
                  sigmoid((y - 30.0 + 5.0)));
  });
  return 0;
}
```

---

# Example pt. 2: Script (`test.lmp`)

```lammps
# Simulation setup goes here

# Add the external controller w/ maxdelay of 50 ms
fix myarbfn all arbfn maxdelay 50.0

# Run the simulation
reset_timestep 0
timestep 0.01
run 4000000
```

---

# Example pt. 3: Compile and run

```sh
# Compile
mpicxx -O3 -std=c++11 -g \
    -o wall_effect.out wall_effect.cpp

# Run
mpirun -n 1 wall_effect.out : \
    -n 3 lmp -mpicolor 123 -in test.lmp
```

---

# Special case: Static force fields

- If we want to impose constant force fields, we can be much
    more efficient
- We provide `fix arbfn/ffield` for this purpose
- Gets a finite set of "nodes" in space w/ `dfx`, `dfy`, and
    `dfz` data from the MPI controller via the `"gridRequest"`
    special case request
- Trilinearly interpolates in $x/y/z$ space for the force deltas
- This runs **as fast as normal LAMMPS**, while still allowing
    for arbitrary force fields

---

# Gallery

---

# Code tour

[Codebase on GitHub](https://github.com/jorbDehmel/lammps_arb_fn)
