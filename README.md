# Introduction

The `GFM_DE_rule` package computes the gravitational potential of a body, 
the shape of which is given by a surface spherical harmonic expansion and
the density of which is constant. The evaluation points can be located inside
the body, on its surface or outside the body.

The package builds on the spatial-domain gravity forward modelling technique
developed by [Fukushima (2017)](https://doi.org/10.3847/1538-3881/aa88b8),
which evaluates the Newton integral in the spatial domain via the double
exponential rule. We used the Fortran routines from that paper, while several
modifications were necessary to fit the goals of the package.

The package is written in Fortran and MATLAB. The Fortran version supports
OpenMP parallelization and can be compiled in double and quadruple precision.
The Fortran code is generally faster than the MATLAB code.

The accuracy of the output potential can reach ~14 correct digits in
double precision and ~30 digits in quadruple precision. The accuracy level
is controlled via a relative tolerance error parameter. The higher the required
accuracy, the longer the computation times can be expected.


# Repository structure

* `./data` gathers all data that are necessary to run test computations (see
  below).

  * The `./data/Bennu_Shape_SHCs_to15.txt` file contains spherical harmonic
    coefficients defining the shape of the 101955 Bennu asteroid up to degree
    `15`. For details on the structure of the file, see, e.g., the
    `./src/Fortran/danm_bnm.f95` subroutine. The development of these
    coefficients is described in the documentation (see the `./docs` folder).

  * `./data/Bennu_Computing_points.txt` specifies three points, at which the
    gravitational potential is computed when running the test computations.

* `./docs` contains a preprint of [Bucha and Sanso
  (2021)](https://doi.org/10.1007/s00190-021-01493-w), in which the GFM_DE_rule
  package was published, and a preprint of [Fukushima
  (2017)](https://doi.org/10.3847/1538-3881/aa88b8), in which the core
  technique and codes of GFM_DE_rule were published.

* `./src` is the directory with the source codes.  The Fortran code is stored
  in `./src/Fortran` and the MATLAB code in `./src/MATLAB`.

  * `./src/Fortran`: The files starting with `d` denote routines that run in
    double precision, while the `q` prefix stands for their quadruple
    counterpart.  Each routine contains a detailed documentation.

  * `./src/MATLAB`: all files have the `d` prefix, given that MATLAB does not
    support quadruple precision.

* `./Test_run.sh` is a very simple Linux bash script to compile and run a test
  computation using the Fortran version of `GFM_DE_rule`.


# Test computation

This section provides a brief guide on how to run the test computations, both
with the Fortran and MATLAB code.

## Fortran code

### Linux

Execute the bash script `./Test_run.sh`. It compiles the test program
`./src/Fortran/dTest_run.f95` (double precision) with GNU Fortran and then
executes the compiled binary.

To run the test computation in quadruple precision
(`./src/Fortran/qTest_run.f95`), set the variable `precision` in
`./Test_run.sh` to `q` and execute the bash script `./Test_run.sh`.

By default, `./Test_run.sh` compiles the code without the flag enabling
OpenMP. The Fortran subroutines, however, support the OpenMP parallelization by
default, so you may want to manually add the `-fopenmp` flag in `./Test_run.sh`
to enable OpenMP.

### Other operating system

For the time being, you have to compile the package manually by yourself.


## MATLAB code

Go to the `./src/MATLAB` folder and execute the `dTest_run.m` script.

To allow parallel execution, open a parpool session in MATLAB (command
`parpool`) before running the script.


# Contributing

Contributions of any kind are welcome!


# Contact

Feel free to contact the author, Blazej Bucha, at blazej.bucha@stuba.sk.


# Citing

The `GFM_DE_rule` package was published by

 * Bucha B, Sanso F (2021). Gravitational field modelling near irregularly 
   shaped bodies using spherical harmonics: a case study for the asteroid
   (101955) Bennu. doi:
   [https://doi.org/10.1007/s00190-021-01493-w](https://doi.org/10.1007/s00190-021-01493-w).

The package is based on the gravity forward modelling method and the Fortran
routines developed by

 * Fukushima T (2017). Precise and fast computation of the gravitational field 
   of a general finite body and its application to the gravitational study of 
   asteroid Eros. The Astronomical Journal 154(145):15pp, doi:
   [https://doi.org/10.3847/1538-3881/aa88b8](https://doi.org/10.3847/1538-3881/aa88b8).

Please consider acknowledging these references when using the package.

