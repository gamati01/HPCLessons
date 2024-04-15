#
rm -rf *.o mm.x *.mod
#
COMP=nvfortran
OPT="-Minfo=acc -O2 -acc"
echo "compiling with " $COMP $OPT 
for i in 0 1 2; do
    echo "compilin mm.$i.F90"
    $COMP $OPT mod_tools.F90 -c
    $COMP $OPT mm.$i.F90 -c
    $COMP $OPT mod_tools.o mm.$i.o -o mm.$i.x
done
echo "That's all folks!!!"

