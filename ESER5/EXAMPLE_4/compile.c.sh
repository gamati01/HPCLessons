#
rm -rf *.o mm.x 
#
COMP=nvcc
#OPT=" -O2 -DSINGLEPRECISION"
OPT=" -O2 "
echo "compiling with " $COMP $OPT 
$COMP $OPT mm.cu -o mm.x
#
OPT=" -O2 -DVALIDATION"
echo "compiling with " $COMP $OPT 
$COMP $OPT mm.cu -o mm.validation.x
echo "That's all folks!!!"

