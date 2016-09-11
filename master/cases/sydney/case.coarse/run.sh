#!/bin/sh
rm -r processor*
decomposePar
mpirun -np 8 lowMachDieselFoam -parallel > log &
