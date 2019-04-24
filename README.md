# GROMACS (Fluctuating Hydrodynamics / Molecular Dynamics two-way coupling)

# Installation
```bash
mkdir build
cd build
cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON
(mb-pro$ cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++)
make
make check # optional
sudo make install
/usr/local/gromacs/bin/gmx
(mb-pro$ /usr/local/gromacs/bin/gmx)
```

# Compile and run

```
make # or
make -j 4 # 4 processor cores
sudo make install
/usr/local/gromacs/bin/gmx grompp
/usr/local/gromacs/bin/gmx mdrun -ntomp 1 # or
/usr/local/gromacs/bin/gmx mdrun -nt 1 -v # -nt 1 for single thread and -v for outputing steps
vmd (New molecule -> confout.gro, traj_comp.xtc)
```
# Open project in Eclipse

Prepare for release




1. Start Eclipse in _Terminal_ ```./eclipse```
1. Build the _Debug_ version instead of _Release_ using ```ccmake gromacs_fhmd``` (Set the option to ```Debug```)
OR
```
cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_BUILD_TYPE=Debug
```

2. Import Existing Makefile Project
3. Configure Build command, _Project > properties > C/C++ Build > Build command_ ```make -j 4```. Set Build directory to ```/home/aware/Desktop/gromacs_fhmd_debug```. Set working directory in tab Arguments to  ```/home/aware/Desktop``` (the path where all the *gmx* files are).
4. Debug configuration, set _C/C++ Application_ path to ```/home/aware/Desktop/gromacs_fhmd_debug/bin/gmx``` with arguments ```mdrun -nt 1```. Set working directory in tab Arguments to  ```/home/aware/Desktop``` (the path where all the *gmx* files are).

# Besides

How to downgrade to gdb 8.0.1

Unlink current gdb: 
```brew unlink gdb
```

Install gdb 8.0.1: 
```brew install https://raw.githubusercontent.com/Homebrew/homebrew-core/9ec9fb27a33698fc7636afce5c1c16787e9ce3f3/Formula/gdb.rb
```
Optional: avoid upgrade gdb with 
```brew pin gdb```

## Signing `gdb`

# Credits
Ivan Korotkin
