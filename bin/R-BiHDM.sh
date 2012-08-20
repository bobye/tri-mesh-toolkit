#!/bin/sh

base=$(basename $1 .off)

./MeshTK -i $base -l
mkdir $base.ev
./eigen_solver -matrix $base.stiff -mass_matrix $base.mass -eps_nev 50 -shift_val 0.01 -output_file $base.ev -eps_monitor
./MeshTK -i $base -mgf
./eigen_solver -matrix $base.fbihdmat -eps_nev 23
