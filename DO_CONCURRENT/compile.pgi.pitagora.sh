#
module purge
module load nvhpc/24.9
module li
#
export OPT=-O3 
rm -rf *.o mm.x *.mod
#
echo "compiling with" $OPT
echo "serial" 
echo "---------------------------------------------------------" 
nvfortran -Minfo=all                   $OPT mod_tools.F90 -c
nvfortran -Minfo=all                   $OPT -DVALIDATION mm.F90 -c
nvfortran -Minfo=all                   $OPT -DVALIDATION mod_tools.o mm.o -o mm.serial.x
#
echo "multicore" 
echo "---------------------------------------------------------" 
nvfortran -stdpar=multicore -mp -Minfo=all $OPT mod_tools.F90 -c
nvfortran -stdpar=multicore -mp -Minfo=all $OPT -DVALIDATION mm.F90 -c
nvfortran -stdpar=multicore -mp -Minfo=all $OPT -DVALIDATION mod_tools.o mm.o -o mm.multi.x
#
echo "gpu" 
echo "---------------------------------------------------------" 
nvfortran -stdpar=gpu       -mp -Minfo=all $OPT mod_tools.F90 -c
nvfortran -stdpar=gpu       -mp -Minfo=all $OPT -DVALIDATION mm.F90 -c
nvfortran -stdpar=gpu       -mp -Minfo=all $OPT -DVALIDATION mod_tools.o mm.o -o mm.gpu.x
#
echo "That's all folks!!!"
