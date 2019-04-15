# GROMACS (Fluctuating Hydrodynamics / Molecular Dynamics two-way coupling)

# Installation
```
mkdir build
cd build
cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON
make
make check
sudo make install
gmx
```

# Compile and run

```
make # or
make -j 4 # 4 processor cores
sudo make install
/usr/local/gromacs/bin/gmx grompp
/usr/local/gromacs/bin/gmx mdrun -ntomp 1 # or
/usr/local/gromacs/bin/gmx mdrun -nt 1 -v # -nt 1 for single process and -v for outputing steps
vmd (New molecule -> confout.gro, traj_comp.xtc)
```
# Open project in Eclipse

1. Start Eclipse in _Terminal_ ```./eclipse```
1. Build the _Debug_ version instead of _Release_ using ```ccmake gromacs_fhmd``` (Set the option to ```Debug```)
2. Import Existing Makefile Project
3. Configure Build command, _Project > properties > C/C++ Build > Build command_ ```make -j 4```. Set Build directory to ```/home/aware/Desktop/gromacs_fhmd_debug```
4. Debug configuration, set _C/C++ Application_ path to ```/home/aware/Desktop/gromacs_fhmd_debug/bin/gmx``` with arguments ```mdrun -nt 1```.

# Credits
Ivan Korotkin
