#!/bin/bash

make clean
make
./bin/lj_gas.out
python analyze_results/read_energies.py
export LD_PRELOAD=/lib64/libsqlite3.so.0
vmd -e show_traj.tcl
