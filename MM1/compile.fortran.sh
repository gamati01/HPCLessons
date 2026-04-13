#
export OPT=$1
#
rm -rf *.o mm.x *.mod
#
# nvfortran (NVIDIA)
COMP=nvfortran
#
# ifx (intel)
COMP=ifx
#
# gfortran (GNU)
COMP=gfortran
#
echo "compiling with " $COMP $OPT 
#
$COMP $OPT mod_tools.F90 -c
$COMP $OPT mm.F90 -c
$COMP $OPT mod_tools.o mm.o -o mm.x
#
echo "That's all folks!!!"
#

