#!/bin/bash
#PBS -l nodes=1:ppn=24
#PBS -l walltime=02:00:00
#PBS -N DIM-OPT

cd "$PBS_O_WORKDIR"
cp ../INCAR_OPT INCAR
cp ../POTCAR .
cp POSCAR POSCAR.F
rm pcell.bin
rm WAVECAR
python /home/lzhu/apps/writekp.py 0.05
mpirun -n 24 /home/lzhu/apps/vasp.5.4.1/bin/vasp_std > log
python ../ovjob.py
