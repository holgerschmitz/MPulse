#!/bin/bash

# processes used
NP=2

# executable path
EXEC="./../../../build/bin/mpulse2d"

# Run the mpirun command
mpirun -np $NP $EXEC

