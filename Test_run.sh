#!/usr/bin/env bash


# This is a simple Linux bash script to compile the Fortran source codes of the
# "GFM_DE_rule" package and to execute the compiled binary.  "gfortran" must be
# installed on your machine.  To run the test computation in MATLAB, execute
# the "./src/MATLAB/dTest_run.m" script in MATLAB.


# Set the "precision" variable to "d" for double precision computation or to
# "q" for quadruple precision computation. The rest of this script may not need
# to be modified.
precision="d"


# Exit when any command fails
set -e


echo "Starting the execution of the \"Test_run.sh\" script..."


# Specify the Fortran files to be compiled
echo "Reading the names of the Fortran files to be compiled"
f_src="./src/Fortran/vartypes.f90
       ./src/Fortran/constants.f90
       ./src/Fortran/${precision}anm_bnm.f90
       ./src/Fortran/${precision}log1p.f90
       ./src/Fortran/${precision}qde1.f90
       ./src/Fortran/${precision}qde2.f90
       ./src/Fortran/${precision}Test_run.f90"


# Specify the name of the executable file to be compiled. The file will be
# saved in the "./bin" folder, which will be automatically created if it does
# not exist.
echo "Reading the name of the executable file..."
f_bin=${precision}Test_run


# Specify compiler flags
echo "Reading the compiler flags..."
compiler_flags="-ffast-math -fexpensive-optimizations -O3"


# Create the "./bin" folder to store the executable file
echo "Creating the ./bin folder..."
mkdir -p bin


# Compile using GNU Fortran
echo "Compiling the package..."
gfortran $compiler_flags -o bin/$f_bin $f_src


# Remove the created "*.mod" files
echo "Removing the created *.mod files..."
rm vartypes.mod constants.mod


# Change the working directory to "./bin" (necessary to allow the executable
# file to find input data based on the relative path)
echo "Changing working directory to ./bin..."
cd bin


# Execute the compiled binary file
echo "Executing the compiled binary file (single-threaded computation by\
 default)..."
if [ "$precision" = "q" ]; then
    echo "NOTE: The default relative accuracy for quadruple precision\
 computations is set to 1e-30. The computation may therefore take several\
 minutes, especially for the default single-threaded computation..."
fi

echo "====================================="
echo ""
echo ""
./$f_bin
echo ""
echo ""
echo "====================================="
echo "End of the Test_run.sh script"

