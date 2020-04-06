#!/bin/sh

# This sh file is to help with the compilation of hoomd!

# First off, you'll need these dependencies:
# 1. Python >= 3.5
# 2. NumPy >= 1.7 (once you have python, just $pip install numpy)
# 3. CMake >=2.8.10.1
# 4. C++11 compiler (already exists on a stock mac)

# Set the paths where you'd like hoomd and its python interface to exist
buildPath="/Users/kolbt"                        # EDIT THIS TO DESIRED BUILD PATH!!!
hoomdPyPath="/Users/kolbt/hoomd-blue/build"     # EDIT THIS TO DESIRED PYTHON MODULE PATH!!!

# Go to the directory you want to put hoomd
cd "${buildPath}"
# Clone the most recent hoomd-blue build from github
git clone --recursive https://github.com/glotzerlab/hoomd-blue

# Then you want to make a build directory within hoomd-blue for compilation
cd hoomd-blue
mkdir build
cd build
# There are other command line options for this next line... see documentation
cmake ../ -DCMAKE_INSTALL_PREFIX="${hoomdPyPath}"
make -j4
make install

exit 0
