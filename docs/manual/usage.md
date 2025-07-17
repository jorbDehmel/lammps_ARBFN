
# Using ARBFN

This document details how to write and use controllers via the
ARBFN package.

**Note:** As of v0.3.1, only one controller can exist within a
simulation. Using multiple controllers will result in undefined
behaviour (the MPI communications will almost certainly
collide).

## `fix arbfn`

The LAMMPS side of the fix just sets up the connection to the
controller with some basic protocol settings. The following code
snippet shows a few valid lines.

```lammps
# The base command: This links atom type 1 to the controller,
# sending dipole information
fix name_1 1 arbfn dipole

# This links all atoms to the controller with a timeout of 50ms
fix name_2 all arbfn maxdelay 50.0

# This links all atoms, updating only every 100 timesteps
fix name_3 all arbfn every 100

# All atoms, updating every 100 steps, timeout of 50ms, dipole
fix name_4 all arbfn maxdelay 50.0 every 100 dipole
```

## `fix arbfn/ffield`

The LAMMPS side of the fix just sets up the connection to the
controller with some basic protocol settings. The following code
snippet shows a few valid lines.

```lammps
# 101 by 201 by 301 POINTS in space
fix name_1 all arbfn/ffield 100 200 300

# 1 by 2 by 3 BINS (2 by 3 by 4 points), updating every 25 steps
fix n2 1 arbfn/ffield 1 2 3 every 25

# You can also send the controller dipole information, but not
# change it. The grid will update every 10 steps
fix n3 all 20 20 1 dipole every 10
```

## Special Case: Controllers in `python 3`

**This is the easiest language to implement controllers in.**
As part of this software, we have provided
[`controller.py`](../../tests/controller.py) which exposes a
simplified wrapper `Python 3`.

```py
# A shebang would be ok too

import math
from typing import Tuple
from controller import ffield_controller


MAG: float = 0.01
WELL: float = 20.0


def main() -> None:
    def apply_force(pos) -> Tuple[float, float, float]:
        '''
        Given a 3-tuple pos containing the position in
        simulation space, return a 3-tuple with the x/y/z forces
        for this point.
        '''
        forces = [0.0] * 3
        forces[0] = forces[1] = forces[2] = 0.0
        if pos[0] != 0.0 or pos[1] != 0.0:
            dist: float = math.sqrt(pos[0] ** 2 + pos[1] ** 2)
            force: float = MAG * (WELL - dist) / WELL
            forces[0] = force * pos[0] / dist
            forces[1] = force * pos[1] / dist
        return forces[0], forces[1], forces[2]

    ffield_controller(apply_force)


# Note: The nuances of MPI in python mean that
# "if __name__ == '__main__'" doesn't work. You must manually
# call any main fn like this.

main()

```

Assuming that the above file is called `py_controller.py`, that
there is a properly formatted `input_script.lmp` using
`fix arbfn`, and that `lmp` is a valid command with ARBFN
installed, the following command would run a simulation with 1
controller and 3 worker threads.

```bash
mpirun -n 1 \
    python3 ./py_controller.py \
    : -n 3 \
    lmp -mpicolor 123 -in input_script.lmp
```

More information can be found in PyDoc format
[in the module](../../tests/controller.py). Wrapper functions
for `fix arbfn` (both "independent" and "dependent") are
documented therein.

## Special Case: Controllers in `C++`

**This is the second easiest language for controllers.**
As part of this software, we have provided
[`controller.hpp`](../../tests/controller.hpp) which exposes a
simplified wrapper for `C++`. This is header-only, so no
additional objects need to be included at compile-time.

The following excerpt demonstrates how to write a controller
for `fix arbfn/ffield` in `C++`.

```cpp
#include "lammps_ARBFN/tests/controller.hpp"
#include <boost/json/src.hpp>
#include <cmath>

const static double rad = 20.0;
const static double mag = 0.1;

int main() {
  // Apply outward force within a small circle in the
  // middle
  ffield_controller([](const boost::json::value &atoms,
                        const bool &first, const double pos[3],
                        double forces[3]) {
    // This lambda is called for every point in the grid upon
    // simulation start. The initial atoms are in `atoms`,
    // `first` is true only the first time the lambda is called,
    // `pos` is a 3-tuple containing the x/y/z position in the
    // simulation space, and `forces` is where we WRITE the
    // information we want to put at that point in space. The
    // lambda does not return a value.
    forces[0] = forces[1] = forces[2] = 0.0;
    if (pos[0] * pos[0] + pos[1] * pos[1] <= rad * rad &&
        (pos[0] != 0.0 || pos[1] != 0.0)) {
      const double dist =
          std::sqrt(pos[0] * pos[0] + pos[1] * pos[1]);
      forces[0] = mag * pos[0] / dist;
      forces[1] = mag * pos[1] / dist;
    }
  });
  return 0;
}
```

Any `C++` controllers should be compiled with `mpicxx` rather
than a non-MPI compiler. More information can be found in
Doxygen format [in the header file](../../tests/controller.hpp).
Wrapper functions for `fix arbfn` (both "independent" and
"dependent") are documented therein.

```bash
# If the above file was controller.cpp, this would compile it to
# controller.out
mpicxx -std=c++11 -o controller.out controller.cpp
```

Assuming that there is a properly formatted `input_script.lmp`
using `fix arbfn`, and that `lmp` is a valid command with ARBFN
installed, the following command would run a simulation with 1
controller and 3 worker threads.

```bash
mpirun -n 1 \
    ./controller.out \
    : -n 3 \
    lmp -mpicolor 123 -in input_script.lmp
```

## Controllers in Any Other Language

If you want to write a controller in any non-listed language,
you will have to handle the MPI communication protocol. This is
easiest to do by translating the existing
[`python`](../../tests/controller.py) or
[`C++`](../../tests/controller.hpp) protocols into your desired
language, as well as by studying the protocol documentation for
[`fix arbfn`](../joss/arbfn_protocol.png),
[`fix arbfn/ffield`](../joss/ffield_protocol.png), and
[`fix arbfn/ffield every n`](../joss/ffield_every_protocol.png).

**Your controller language must be able to handle MPI** and
should probably be able to handle JSON. Any language that meets
these criteria (and is sufficiently nontrivial) can be used to
create a controller.
