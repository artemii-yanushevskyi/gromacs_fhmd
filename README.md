# GROMACS (Fluctuating Hydrodynamics / Molecular Dynamics two-way coupling)

Modified by Ivan Korotkin.

# Installation

mkdir build
cd build
cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON
make
make check
sudo make install
gmx
