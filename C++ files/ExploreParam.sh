# Basic while loop
#           NLOCI   NPLOIDY     NINIT[0]    NINIT[1]    NLOCAL      SC_MAJOR    SC_LOCAL    NGEN    NREP    RECOMB  K

#!/bin/bash

for nin1 in 10
do
for dist in 1 2 3 4 5
do
for rec in 0.0 0.005 0.01 0.015 0.02 0.025 0.03 0.035 0.04 0.045 0.05 0.055 0.06 0.065 0.07 0.075 0.08 0.085 0.09 0.095 0.1 0.2 0.3 0.4 0.5
do
./test.out  50      2           200        20          2       1.2         1.0        75      200     0.01      200
done
done
done

./test.out  50      2           200        40          20       0.2         -0.05        200      200     0.01      250
