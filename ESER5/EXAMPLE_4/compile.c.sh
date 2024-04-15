#
rm -rf *.o mm.x 
#
COMP=nvcc
#OPT=" -O2 -DSINGLEPRECISION"
OPT=" -O2 "
echo "compiling with " $COMP $OPT 
$COMP $OPT mm.cu -o mm.x
echo "That's all folks!!!"

