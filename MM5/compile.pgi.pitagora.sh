#
#
module purge
module load nvhpc/24.9
module li
#
export COMP=nvfortran
export OPT="-mp=gpu -Minfo=all "
rm -rf *.o mm.x *.mod
#
echo "compiling with" $COMP $OPT $i
$COMP $OPT mod_tools.F90 -c
$COMP $OPT mm.F90 -c
$COMP $OPT mod_tools.o mm.o -o mm.fast.x
#
$COMP $OPT -DVALIDATION mod_tools.F90 -c
$COMP $OPT -DVALIDATION mm.F90 -c
$COMP $OPT -DVALIDATION mod_tools.o mm.o -o mm.check.x
echo "That's all folks!!!"

