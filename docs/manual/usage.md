
# Using ARBFN

This document details how to write and use controllers via the
ARBFN package.

## `fix arbfn`

**TODO: WRITE THIS**

## `fix arbfn/ffield`

**TODO: WRITE THIS**

## Special Case: Controllers in `python 3`

**This is the easiest language to implement controllers in.**
As part of this software, we have provided
[`controller.py`](../../tests/controller.py) which exposes a
simplified wrapper `Python 3`.

```py
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


main()

```

More information can be found in PyDoc format
[in the module](../../tests/controller.py). Wrapper functions
for `fix arbfn` (both "independent" and "dependent") are
documented therein.

## Special case: Controllers in `C++`

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

## Controllers in Any Other Language

If you want to write a controller in any non-listed language,
you will have to handle the MPI communication protocol. This is
easiest to do by translating the existing `Python` or `C++`
protocols into your desired language, as well as by studying the
protocol documentation for
[`fix arbfn`](../joss/arbfn_protocol.png),
[`fix arbfn/ffield`](../joss/ffield_protocol.png), and
[`fix arbfn/ffield every n`](../joss/ffield_every_protocol.png).

**Your controller language must be able to handle MPI** and
should probably be able to handle JSON. Any language that meets
these criteria (and is sufficiently nontrivial) can be used to
create a controller.
