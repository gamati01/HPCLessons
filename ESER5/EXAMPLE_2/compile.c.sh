#
rm -rf *.o mm.x *.mod
#
COMP=nvc
OPT="-Minfo=acc -O2 -acc"
echo "compiling with " $COMP $OPT 
$COMP $OPT mm.c -c
$COMP $OPT mm.o -o mm.x
echo "That's all folks!!!"

