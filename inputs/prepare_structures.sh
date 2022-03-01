#!/bin/sh
#BSUB -q bowman
#BSUB -gpu "num=1:gmodel=TeslaP100_PCIE_16GB"
#BSUB -n 14
#BSUB -R span[ptile=14]
#BSUB -J equil
#BSUB -o equil-%J.log
#BSUB -e equil-%J.log

source /project/bowmore/mizimmer/installations/gromacs_2020.1_p100_build

d='/project/bowmanlab/jjmiller/Simulations/PlasmepsinV/SetupFiles'
prot='delNep-prot-masses.pdb'

echo 1 | gmx pdb2gmx -f $prot -o b4editconf.gro -p topol.top -water tip3p -ignh -vsite hydrogens -ff amber03
#echo 1 | gmx pdb2gmx -f $prot -o b4editconf.gro -p topol.top -water tip3p -ignh
gmx editconf -c -f b4editconf.gro -o b4solv.gro -bt dodecahedron -d 1.0
gmx solvate -cp b4solv.gro -cs -o b4ions.gro -p topol.top

mv topol.top topol-noions.top
gmx grompp -f $d/em.mdp -c b4ions.gro -p topol-noions.top -o ions
echo 'SOL' | gmx genion -s ions.tpr -o b4em.gro -p topol-noions.top -neutral -conc 0.100
mv topol-noions.top topol.top

gmx grompp -f $d/em.mdp -c b4em.gro -p topol.top -o em
gmx mdrun -s em.tpr -o em -c b4equil -v -nt 8 -ntmpi 1 -ntomp 8

gmx grompp -f $d/pr.mdp -c b4equil.gro -p topol.top -o md -maxwarn 1 -r b4equil.gro
gmx mdrun -s md -o md -c after_equil -v -nt 8 -ntmpi 1 -ntomp 8 -x after_equil -pin on
# align trajectory and starting structure (optional and make sure this doesn't do anything weird with periodic boundaries)
echo '10 0' | gmx trjconv -f after_equil.xtc -o after_equil_aligned.xtc -s md.tpr -pbc mol -ur compact -center
echo '10 0' | gmx trjconv -f after_equil.gro -o start.gro -s md.tpr -pbc mol -ur compact -center

