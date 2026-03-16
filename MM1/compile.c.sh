#
rm -rf *.o mm.x *.mod
#
# nvc compiler (nvidia)
COMP=nvc
#OPT=
#
# icx compiler (intel)
COMP=icx
#OPT=
#
# gcc compiler (GNU)
COMP=gcc
OPT=
#
echo "compiling with " $COMP $OPT 
#
$COMP $OPT mm.c -c
$COMP $OPT mm.o -o mm.x
#
echo "That's all folks!!!"
#
