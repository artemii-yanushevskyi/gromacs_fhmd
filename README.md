# GROMACS (Fluctuating Hydrodynamics / Molecular Dynamics two-way coupling)

Modified by Ivan Korotkin.

# Energy Minimisation

	gmx grompp -f minim.mdp -c conf.gro -p topol.top -o em.tpr -maxwarn 1
	gmx mdrun -v -deffnm em
	
	(original)gmx grompp -f grompp.mdp -c em.gro -o afterminim_no_fh.tpr
	(original)gmx mdrun -nt 1 -s afterminim_no_fh.tpr -c aftermimim_no_fh.gro 

	gmx grompp -f grompp.mdp -c aftermimim_no_fh.gro -o afterminim_fh.tpr
	gmx mdrun -nt 1 -s afterminim_fh.tpr -o afterminim_traj.trr -c afterminim_confout.gro -v
