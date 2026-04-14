#
module load nvhpc/23.1
module li
#
rm -rf *.o mm.x *.mod
#
# gfortran
COMP=gfortran
OPT="-fopenmp -O3" 
#
# intel
COMP=ifx
OPT="-qopenmp -O3" 
#
# nvidia
COMP=nvfortran
OPT="-Minfo=mp -O3 -mp" 
#
echo "compiling with " $COMP $OPT 
$COMP $OPT mod_tools.F90 -c
$COMP $OPT mm.F90 -c
$COMP $OPT mod_tools.o mm.o -o mm.x
echo "That's all folks!!!"

