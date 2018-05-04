#!/bin/bash
set -e
make -f Makefile
module load buildenv-intel/2015-1
usereservation tddc78-2018-$1

for i in {1..4}; do
	cores=$( echo "2^$i" | bc )
	echo "Cores: $cores \n" >> "result$i.txt"
	salloc -N1 -n$cores sh test.sh $2 $i
done
