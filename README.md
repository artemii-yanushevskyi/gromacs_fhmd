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

# Compile and run


make
sudo make install
/usr/local/gromacs/bin/gmx grompp
/usr/local/gromacs/bin/gmx mdrun -ntomp 1
vmd (New molecule -> confout.gro, traj_comp.xtc)
