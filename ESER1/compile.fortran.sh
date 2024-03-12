#
rm -rf *.o mm.x *.mod
#
# nvfortran (NVIDIA)
COMP=nvfortran
OPT=
#
# ifx (intel)
COMP=ifx
OPT=
#
# gfortran (GNU)
COMP=gfortran
OPT=-O3
#
echo "compiling with " $COMP $OPT 
#
$COMP $OPT mod_tools.F90 -c
$COMP $OPT mm.F90 -c
$COMP $OPT mod_tools.o mm.o -o mm.x
#
echo "That's all folks!!!"
#

