#!/bin/bash

make clean
make
./bin/lj_gas.out
python analyze_results/read_energies.py
vmd -e show_traj.tcl
