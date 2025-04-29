#
rm -rf *.o mm.x *.mod
#
# version
#
VER=4
#
# nvfortran (NVIDIA)
COMP=nvfortran
OPT=-O2
#
# ifx (intel)
#COMP=ifx
#OPT=
#
# gfortran (GNU)
COMP=gfortran
OPT=-O2
#
echo "compiling, mm.$VER.F90 with " $COMP $OPT 
#
$COMP $OPT mod_tools.F90 -c
$COMP $OPT mm.$VER.F90 -c
$COMP $OPT mod_tools.o mm.$VER.o -o mm.$VER.x
#
echo "That's all folks!!!"
#

