#!/bin/bash

# processes used
NP=4

# executable path
EXEC="./../../../build/bin/mpulse2d"

# Run the mpirun command
mpirun -np $NP $EXEC

