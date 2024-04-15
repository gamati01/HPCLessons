for size in 1024 2048 4096 8192; do
	echo $size > input.dat
	echo $size 
	./mm.x < input.dat > output.$size.log
done
