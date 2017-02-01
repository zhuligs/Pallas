#!/bin/bash
#PBS -l nodes=1:ppn=24
#PBS -l walltime=02:00:00
#PBS -N DIM-TS

cd "$PBS_O_WORKDIR"
rm WAVECAR
rm pcell.bin
cp PRESAD.vasp POSCAR.F
python -u ../dvjob.py > dvjob.out