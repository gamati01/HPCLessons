#
#!/bin/bash
#
THREADS=$1
#
if [ -z "$THREADS" ]; then
    echo "Threads not set "
    echo "Usage: $0 <threads>"
    exit 1
fi
#
#
for size in 2048 4096 8192; do
	echo "Scaling up to $1 threads for size $size"
	echo $size > input.dat
#
	for th in `seq 1 $1`; do
		export OMP_NUM_THREADS=$th
		echo "threads = "$th
		./mm.x < input.dat > out.$size.$th.dat
	done	
done
