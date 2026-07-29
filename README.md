# MPulse

MPulse is a parallel C++ finite-difference time-domain (FDTD) solver for the full Maxwell equations in one, two, and three spatial dimensions. It is intended for simulations of electromagnetic pulses in linear, nonlinear, dispersive, and ionising media.

The active code supports plain Yee-style FDTD, instantaneous and time-averaged Kerr nonlinearities, CPML absorbing boundaries, externally injected fields, damped plasma currents, spatially varying plasma densities, and experimental multiphoton and field-ionisation routines. A three-pole piecewise-linear recursive-convolution (PLRC) implementation is also present in the source tree, but it is currently disabled and is not part of the default build.

> **Repository status**
>
> This is research code rather than a packaged application. The supplied archive is not self-contained: the top-level `Makefile` expects a populated `huerto/` source tree, while the archive contains only an empty `huerto/` directory. The build file also contains machine-specific include and library paths. See [Building](#building) and [Known limitations](#known-limitations).

## Features

- Full-vector Maxwell solver on a staggered Yee grid.
- Separate 1D, 2D, and 3D executables selected at compile time.
- MPI Cartesian domain decomposition and halo exchange.
- Plain FDTD through `FDTD_Plain`.
- Instantaneous Kerr response through `FDTD_Kerr`:

  \[
  \mathbf D = \varepsilon_0\left(\varepsilon_r + \chi |\mathbf E|^2\right)\mathbf E.
  \]

- Time-averaged Kerr response through `FDTD_KerrAverage`.
- Drude-like damped plasma current coupled to a density field.
- Spatially varying electron and neutral density fields initialised from setup-file expressions.
- Experimental multiphoton ionisation and strong-field/ADK-like ionisation tasks.
- Convolutional perfectly matched layer (`CPMLBorder`).
- Plane-wave, plane-Gaussian, and Gaussian-beam sources.
- Three-dimensional ultrashort-pulse source implementation.
- HDF5 field and slice diagnostics.
- Automatic generation of `information.tex` and `references.bib` for numerical-method references used in a run.
- Doxygen configuration for API documentation.

## Repository layout

```text
MPulse/
├── src/             Main MPulse solvers and additional physics
├── huerto/          Required Huerto framework sources; empty in this archive
├── examples/        1D, 2D, and 3D setup-file examples
├── tools/           Legacy HDF5/conversion utilities
├── test/            Legacy tests and test data
├── src_old/         Older implementations retained for reference
├── src_next/        Work-in-progress source implementations
├── Makefile         Builds mpulse1d, mpulse2d, and mpulse3d
├── Doxyfile         Doxygen configuration
└── LICENSE          GNU GPL version 3
```

The production build uses `src/` and selected files under `huerto/`. Files under `src_old/` and `src_next/` are not included by the top-level `Makefile`.

## Requirements

The current Makefile and includes indicate the following dependencies:

- A C++17 compiler.
- An MPI implementation with a C++ compiler wrapper, such as `mpicxx`.
- The Schnek simulation library and headers.
- The Huerto source tree placed under `huerto/` with the directory structure expected by the Makefile.
- HDF5 development headers and libraries.
- FFTW3 development headers and libraries.
- Boost headers; some utilities also require Boost.Program_options and Boost.Format.
- GNU Make.
- Doxygen, optionally, for API documentation.

All dependencies must be ABI-compatible and visible to the selected MPI compiler wrapper.

### The missing `huerto/` tree

The solver includes files such as:

```text
huerto/simulation/simulation_context.hpp
huerto/electromagnetics/fdtd/fdtd_plain.cpp
huerto/electromagnetics/pml/cpml_border.cpp
huerto/diagnostic/field_diagnostic.hpp
```

These files are referenced by both the source and the Makefile, but they are absent from the supplied archive. Restore the Huerto directory from the original checkout, submodule, or source distribution before building.

## Building

The default Makefile uses `mpiCC` and contains site-specific paths under `/home/vol07/...`. Override these variables for your system or edit the Makefile.

From the repository root:

```bash
make \
  CXX=mpicxx \
  LINK=mpicxx \
  INCLUDE="-I/path/to/schnek/include -I/path/to/hdf5/include" \
  LDFLAGS="-L/path/to/schnek/lib -L/path/to/hdf5/lib -L/path/to/fftw/lib"
```

The link libraries in the current Makefile are:

```text
-lhdf5 -lschnek -lfftw3 -lm
```

A successful build creates:

```text
bin/mpulse1d
bin/mpulse2d
bin/mpulse3d
```

To remove generated objects and executables:

```bash
make clean
```

The uploaded archive cannot be build-verified as supplied because its required `huerto/` sources are missing.

## Running a simulation

MPulse has no command-line option for selecting the input file. Each executable opens a file named `mpulse.setup` in the current working directory.

Choose an executable whose dimensionality matches the setup file. For example:

```bash
cd examples/2dimensional/gaussian_beam
mpirun -np 4 ../../../bin/mpulse2d
```

For a serial run through MPI:

```bash
mpirun -np 1 ../../../bin/mpulse2d
```

During a run, the master rank prints the physical time and timestep. Diagnostics write HDF5 files in the working directory. The master rank also writes:

```text
information.tex
references.bib
```

These files record literature references registered by the numerical components used in the simulation.

## Setup-file overview

MPulse uses Schnek's block-oriented setup language. Setup files can define scalar expressions, use mathematical functions and constants, initialise fields from coordinates, and nest source, current, and boundary blocks inside a field solver.

Common built-in constants used by the examples include `pi` and `clight`. The examples use SI-scale lengths, times, electric fields, number densities, electron charge, and electron mass.

### Global parameters

| Parameter | Meaning | Default or requirement |
|---|---|---|
| `Nx`, `Ny`, `Nz` | Number of cells in each compiled dimension | Default array value is `100` |
| `Lx`, `Ly`, `Lz` | Physical domain length in each dimension | Required |
| `tMax` | End time | Required |
| `cflFactor` | Factor used to construct the timestep | `0.99` |
| `ignore_initial_time_stagger` | Skip the initial half magnetic-field step when set to `1` | `0` |

The code computes

```text
dx[d] = L[d] / N[d]
dt    = cflFactor * min(dx) / clight
```

For multidimensional vacuum FDTD, the setup examples choose dimension-dependent factors such as `1/sqrt(2)` in 2D and `1/sqrt(3)` in 3D. Stability remains the user's responsibility when materials, currents, or nonlinearities impose tighter constraints.

### Main blocks

| Block | Purpose | Status |
|---|---|---|
| `EMFields` | Initialises `Ex`, `Ey`, `Ez`, `Bx`, `By`, and `Bz` from expressions | Active; created automatically if omitted |
| `FDTD_Plain` | Linear full-Maxwell FDTD solver | Active |
| `FDTD_Kerr` | Instantaneous isotropic Kerr solver | Active |
| `FDTD_KerrAverage` | Kerr solver using averaged squared-field components | Active |
| `FDTD_PLRC` | Linear three-pole Lorentz/PLRC solver | Source present, disabled |
| `FDTD_PLRC_Nonlinear` | Nonlinear PLRC solver | Source present, disabled |
| `CPMLBorder` | Convolutional PML boundary current | Active |
| `PlaneWaveSource` | Plane-wave injection | Active |
| `PlaneGaussSource` | Plane source with Gaussian pulse envelope | Active |
| `GaussBeamSource` | Gaussian beam injection | Active in 2D/3D builds |
| `ShortPulseInject` | Ultrashort focused-pulse injection | Implemented only in 3D; currently attached only to disabled PLRC solvers |
| `PlasmaDensity` | Initial electron-density field named `Rho` | Active |
| `NeutralDensity` | Initial neutral-density field named `NeutralDensity` | Active |
| `PlasmaCurrent` | Damped plasma current coupled to `Rho` and the electric field | Active as a solver child |
| `MPIIonization` | Multiphoton-ionisation task | Experimental |
| `FieldIonization` | Strong-field/ADK-like ionisation task | Experimental |
| `FieldDiag` | Full-field HDF5 output at a physical-time cadence | Active |
| `GridDiag` | Generic grid HDF5 output | Active |
| `SliceDiag` | Lower-dimensional HDF5 slice at an iteration cadence | Active |

Parameters belonging to Huerto components, including CPML and standard sources, are best inferred from the examples until the missing Huerto documentation/source tree is restored.

## Minimal 1D example

Create `mpulse.setup` in a run directory:

```text
float lambda = 1e-6;
float dx = 0.05*lambda;

Nx = 400;
Lx = Nx*dx;

cflFactor = 1.0;
float dt = cflFactor*dx/clight;
tMax = 3*Lx/clight;

EMFields { }

FDTD_Plain {
  PlaneGaussSource {
    length = 10*lambda;
    originx = -Lx;
    kx = 2*pi/lambda;

    Hx = 0.0;
    Hy = 1.0;
    Hz = 0.0;
  }
}

FieldDiag EyOutput {
  file = "Ey_#t.h5";
  field = "Ey";
  deltaTime = 10*dt;
}
```

Run it with:

```bash
/path/to/MPulse/bin/mpulse1d
```

or, when required by the MPI installation:

```bash
mpirun -np 1 /path/to/MPulse/bin/mpulse1d
```

## Electromagnetic solvers

### `FDTD_Plain`

The plain solver is supplied by Huerto and advances the electric and magnetic fields on a staggered grid. It can contain source, current, and CPML child blocks.

```text
FDTD_Plain {
  PlaneWaveSource { ... }
  PlasmaCurrent { ... }
  CPMLBorder { ... }
}
```

### `FDTD_Kerr`

The instantaneous Kerr solver accepts:

| Parameter | Meaning | Default |
|---|---|---|
| `eps` | Linear relative permittivity | `1.0` |
| `chi` | Cubic nonlinear coefficient used in `(eps + chi*E^2)E` | `0.0` |

At each electric-field update, the solver forms the displacement magnitude and uses a Newton iteration to recover the electric-field magnitude.

```text
FDTD_Kerr {
  eps = 1.0;
  chi = 0.005;
  CPMLBorder { ... }
}
```

### `FDTD_KerrAverage`

The averaged Kerr solver accepts:

| Parameter | Meaning | Default |
|---|---|---|
| `eps` | Linear relative permittivity | `1.0` |
| `chi` | Cubic nonlinear coefficient | `0.0` |
| `T` | Averaging timescale | `-1.0` |
| `E2` | Initial expression for the averaged squared-field components | `0.0` |

It registers diagnostic fields named `E2x`, `E2y`, and `E2z`.

## Plasma model

### Density fields

`PlasmaDensity` creates the field `Rho`, interpreted by `PlasmaCurrent` as a particle number density. Its value can be any setup expression involving coordinates and previously defined variables.

```text
PlasmaDensity {
  float rs = (x-Lx/2)/lambda;
  Rho = 5e26*exp(-rs*rs);
}
```

`NeutralDensity` creates a separate field named `NeutralDensity` and accepts an additional parameter `A`:

```text
NeutralDensity {
  A = 1.0;
  Rho = 5e26;
}
```

Despite the parameter name `Rho` inside the block, the registered neutral field is named `NeutralDensity`.

### Damped plasma current

`PlasmaCurrent` is a child of an electromagnetic solver and registers `PlasmaJx`, `PlasmaJy`, and `PlasmaJz`. The implemented update corresponds to a damped current driven by the electric field and local density:

\[
\frac{\partial \mathbf J}{\partial t} + \gamma\mathbf J
= \frac{Z^2 q^2}{m}\,\rho\mathbf E.
\]

| Parameter | Meaning | Default |
|---|---|---|
| `charge` | Particle charge magnitude | `1.602176634e-19` |
| `mass` | Particle mass | `9.1093837015e-31` |
| `gamma` | Collision/damping frequency | `0.01` |
| `Z` | Charge-state multiplier | `1.0` |

Example:

```text
PlasmaDensity {
  Rho = 5e26;
}

FDTD_Plain {
  PlasmaCurrent {
    charge = 1.602176634e-19;
    mass = 9.1093837015e-31;
    Z = 1.0;
    gamma = 1/(100e-15);
  }
}
```

`PlasmaCurrent` requires a registered `Rho` field. Include a `PlasmaDensity` block before attempting to use it.

## Ionisation models

Ionisation blocks are root-level simulation tasks. They execute after each electromagnetic solver step in the `ionization` task phase and require all of the following fields:

- `Ex`, `Ey`, and `Ez` from `EMFields`.
- `Rho` from `PlasmaDensity`.
- `NeutralDensity` from `NeutralDensity`.

### `MPIIonization`

The multiphoton routine uses the local intensity proxy

```text
I = Ex*Ex + Ey*Ey + Ez*Ez
```

and a power `I^K` in its density update.

| Parameter | Meaning in the current implementation | Default |
|---|---|---|
| `K` | Multiphoton order | `5` |
| `mpa` | Rate coefficient multiplying `I^K` | `1.0` |
| `rat` | Saturation/source density used by the update | `1.0` |
| `Wion` | Used in the plasma-absorption stability warning | `1.0` |

```text
PlasmaDensity { Rho = 0.0; }
NeutralDensity { Rho = 5e26; }

MPIIonization {
  K = 5;
  mpa = 1.0;
  rat = 5e26;
  Wion = 1.0;
}
```

### `FieldIonization`

The field-ionisation routine evaluates an ADK-like rate from the electric-field magnitude.

| Parameter | Meaning in the current source | Default |
|---|---|---|
| `Z` | Ionic charge | `1.0` |
| `Uion` | Ionisation energy, used relative to `13.6` in the rate expression | `1.0` |
| `rat` | Saturation/source density used by the density update | `1.0` |
| `Wion` | Registered but not used by the current rate or density update | `1.0` |
| `neff` | Registered, but the current rate function recomputes a local effective quantum number | `1.0` |

These routines should be treated as experimental. In the current source, both routines retrieve the neutral-density field but update only the electron-density field; they do not deplete `NeutralDensity`. Validate units, coefficients, timestep sensitivity, and particle conservation before using their output quantitatively.

## Sources and boundaries

The example suite demonstrates the active Huerto source blocks:

- `PlaneWaveSource`: monochromatic plane-wave injection with a ramp.
- `PlaneGaussSource`: plane-wave injection with a Gaussian temporal/spatial envelope.
- `GaussBeamSource`: focused Gaussian-beam injection in 2D and 3D.
- `CPMLBorder`: convolutional PML with parameters such as `d`, `sigmaMax`, `kappaMax`, and `aMax`.

Representative CPML configuration:

```text
CPMLBorder {
  d = 10;
  sigmaMax = 8.7;
  kappaMax = 11.0;
  aMax = 0.0;
}
```

Source blocks are nested inside the selected FDTD solver. `GaussBeamSource` is not registered in the 1D build.

## Diagnostics and output

### Full-field diagnostic

```text
FieldDiag EyOutput {
  file = "Ey_#t.h5";
  field = "Ey";
  deltaTime = 10*dt;
}
```

`deltaTime` is a physical-time cadence. The `#t` token in the filename is expanded by the diagnostic framework.

### Slice diagnostic

```text
SliceDiag BzSlice {
  file = "Bz_slice_#t.h5";
  field = "Bz";
  interval = 1000;
  dim = 0;
  pos = 1000;
}
```

`interval` is a timestep cadence. `dim` selects the normal dimension and `pos` selects the grid index of the slice.

Available electromagnetic fields are:

```text
Ex Ey Ez Bx By Bz
```

Additional registered fields can include:

```text
Rho NeutralDensity
PlasmaJx PlasmaJy PlasmaJz
E2x E2y E2z
```

The exact HDF5 layout is controlled by Schnek/Huerto diagnostics.

## Examples

Example setup files are grouped by dimensionality:

```text
examples/1dimensional/
examples/2dimensional/
examples/3dimensional/
```

They cover:

- Plane-wave propagation.
- Explicit plane-wave injection.
- Gaussian pulses.
- Gaussian beams in 2D and 3D.
- Kerr propagation.
- CPML absorption.
- Slice diagnostics.
- A 2D plasma-current case.

Use the matching executable for each directory. Several examples are large enough to require substantial memory, particularly the 3D cases; reduce `Nx`, `Ny`, `Nz`, and the number of diagnostics for initial tests.

## PLRC dispersion status

The files `src/fdtd_plrc.hpp`, `src/fdtd_plrc.cpp`, and `src/fdtd_plrc.t` implement a three-pole Lorentz PLRC model with optional nonlinearity. The intended setup parameters are:

| Parameter family | Meaning |
|---|---|
| `eps` | Infinite-frequency/background relative permittivity |
| `L1eps` ... `L3eps` | Strength of each Lorentz pole |
| `L1delta` ... `L3delta` | Damping parameter for each pole |
| `L1Om2` ... `L3Om2` | Squared resonance-frequency parameter for each pole |
| `chi` | Nonlinear coefficient for the nonlinear PLRC variant |

However, PLRC is disabled in two independent ways:

1. `fdtd_plrc.hpp` is wrapped in `#ifdef DISABLED_CODE`.
2. The include and class registrations in `src/mpulse.cpp` are commented out.

There is also a build-system inconsistency: the top-level `$(wildcard src/*.cpp)` still selects `src/fdtd_plrc.cpp`, even though its class declarations are hidden when `DISABLED_CODE` is not defined. A non-PLRC build should explicitly exclude `src/fdtd_plrc.cpp`; a PLRC development build must instead reconcile the preprocessor guard, registrations, and current Huerto interfaces.

The source itself marks the solver as needing testing. Re-enabling PLRC therefore requires code changes and numerical verification. It should not be assumed to work merely by adding an `FDTD_PLRC` block to a setup file.

Because `ShortPulseInject` is currently registered only as a child of the disabled PLRC solvers, it is not reachable through the active plain or Kerr solver configurations without additional registration changes.

## Documentation

Generate API documentation with:

```bash
doxygen Doxyfile
```

The current Doxygen input is `./src`. The source contains class documentation and citations for the Yee method, PLRC, Kerr/PLRC work, CPML-related components supplied by Huerto, and the ultrashort-pulse source.

## Tests and utilities

The `test/` and `tools/` directories are legacy components and are not integrated into the top-level build.

- `test/Makefile` depends on a user-specific `$(HOME)/src/makeinclude/include.mak` and references an object that is not present in this archive.
- `tools/Makefile` uses the same site-specific include mechanism and copies binaries directly into `$(HOME)/bin/`.
- Utility sources include HDF5-to-text/image conversion, HDF5 joining, slicing, and a 2D FFT tool.

Expect to modernise these Makefiles before using them on a new system.

## Known limitations

- The uploaded archive lacks the required Huerto source tree and cannot be compiled as-is.
- The top-level, test, and tools Makefiles contain machine- or user-specific paths.
- The executable always reads `mpulse.setup` from the current directory; there is no input-file command-line argument.
- PLRC dispersion is present in source form but disabled and explicitly marked as untested; the wildcard source list also needs adjustment for a non-PLRC build.
- `ShortPulseInject` is implemented only for 3D and is currently attached only to disabled PLRC blocks.
- Ionisation routines are experimental, do not currently deplete the neutral-density field, and contain registered parameters that are not used by the implemented update.
- There is no top-level automated test target or continuous-integration configuration in this archive.
- The example problems do not constitute validation for every combination of solver, source, boundary, current, and ionisation model.
- Research-scale 3D grids can require very large amounts of memory and output storage.

## Development notes

When adding a new setup-file block:

1. Implement it as a Schnek block or Huerto simulation entity.
2. Register its class name in `main()` in `src/mpulse.cpp`.
3. Add it as a permitted child of the appropriate parent block.
4. Register shared fields before retrieving them from another component.
5. Decide where it runs in the task or field-solver update sequence.
6. Add a small-dimensional regression example and document its expected result.

Be especially careful with the time staggering between electric fields, magnetic fields, currents, and post-step ionisation. `ignore_initial_time_stagger` is intended for initial data that have already been prepared at the solver's staggered times.

## License

MPulse is distributed under the GNU General Public License, version 3. See [`LICENSE`](LICENSE) for the full terms.
