#
module purge
module load nvhpc/24.9
module li
#
export OPT=-O2 
rm -rf *.o mm.x *.mod
#
echo "compiling with" $OPT
nvfortran -acc -Minfo=all $OPT mod_tools.F90 -c
nvfortran -acc -Minfo=all $OPT -DVALIDATION mm.F90 -c
nvfortran -acc -Minfo=all $OPT -DVALIDATION mod_tools.o mm.o -o mm.x
#
nvfortran -acc -Minfo=all $OPT mod_tools.F90 -c
nvfortran -acc -Minfo=all $OPT mm.F90 -c
nvfortran -acc -Minfo=all $OPT mod_tools.o mm.o -o mm.fast.x
echo "That's all folks!!!"
