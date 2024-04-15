#
rm -rf *.o mm.x *.mod
#
COMP=nvfortran
OPT="-Minfo=acc -O2 -mp=gpu"
echo "compiling with " $COMP $OPT 
$COMP $OPT mod_tools.F90 -c
$COMP $OPT mm.F90 -c
$COMP $OPT mod_tools.o mm.o -o mm.x
echo "That's all folks!!!"

