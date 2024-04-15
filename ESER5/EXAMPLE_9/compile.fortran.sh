#
rm -rf *.o mm.x *.mod
#
COMP=mpif90
OPT="-Minfo -O2 -mp"
echo "compiling with " $COMP $OPT 
$COMP $OPT mm_mpi.F90 -o mm_mpi.x
echo "Launch with: mpirun -np 2 ./mm_mpi.x  2048"

