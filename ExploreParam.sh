# Basic while loop
#           NLOCI   NPLOIDY     NINIT[0]    NINIT[1]    NLOCAL      SC_MAJOR    SC_LOCAL    NGEN    NREP    RECOMB  K

# should look like this ./test.out 20 2 90 10 1 1.1 0.99 100 10 0.5 100
#!/bin/bash
for nin1 in 1 5 10
do
for nloc in 0 1 5
do
for rec in 0.01 0.1 0.5
do
for k in 50 100 200
do
for scm in 1.01 1.05 1.1
do
for scl in 0.99 0.95
do
./test.out  20      2           $k-$nin1    $nin1       $nloc       $scm        $scl        25     10     $rec     $k
done
done
done
done
done
done



