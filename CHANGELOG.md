
# Todo
- Remove `stl` dependencies
- If possible, remove `boost::json` dependency
- Make `ffield_interchange` timed

# Changelog

## `0.3.1` (5/27/2025)
- Made the `every` keyword updatable in `ffield` response
    packets by the controller

## `0.3.0` (3/24/2025)
- **Removed** the `dumponly` flag: There are better ways to do
    this without this package, no need to reinvent the wheel
- Minor `ffield_interchange` API change
- Added `ffield` "refresh every $n$" mode (`every` keyword)
- Updated documentation accordingly

## `0.2.2` (3/3/2025)
- Improved documentation
- Added Doxyfile
- Added check for documentation best practices in `make test`

## `0.2.1` (2/27/2025)
- Added timeout on controller side to allow them to die within
    a set time if the workers fail

## `0.2.0` (2/6/2025)
- Changed API of `ffield` to use node counts instead of bin
    counts, since it's less confusing
- Updated docs accordingly

## `0.1.4` (2/5/2025)
- Added `fix arbfn/ffield` for more efficient constant force
    fields
- Updated docs

## `0.1.3` (2/2/2025)
- Finished presentation in `docs/pres`
- Added some vscode dotfiles
- Removed erroneously added files in `docs/paper`
- Added `docs` rule in root Makefile

## `0.1.2` (1/10/2025)
- Added support for dipole moments (orientation components and
    magnitude $\mu$) via the `dipole` fix argument
- Normalized doxygen in headers

## `0.1.1` (1/9/2025)
- Added some integration testing
- Reorganized repo to unclutter root dir
- Fixed formatting in INSTALL.py
- Removed most output from tests
- Fixed bug causing hang w/ `waiting` packets

## `0.1.0` (1/9/2025)
- Got LAMMPS integration functional
- Documentation improvements
- Added LAMMPS test script
- Added INSTALL.py for quick LAMMPS builds

## `0.0.2` (12/20/2024)
- Got LAMMPS compile functional
- Removed dependency on `boost::mpi` (now just using `boost` for
    JSON)
- Minor documentation improvements

## `0.0.1` (12/16/2024)
- Added minimal functional controller protocol implementations
    in `C++` and `python`
- Wrote out protocol in `README`

## `0.0.0` (12/12/2024)
- Initial commit
