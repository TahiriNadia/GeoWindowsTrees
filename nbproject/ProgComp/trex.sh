#!/bin/bash

declare -i toto; toto=0;
X=0;

mpicc trexParallele.c seq.c phylip.c -o trexPar -lm

rm ~/Parallele/ProgComp/trexPar
cp trexPar ~/Parallele/ProgComp
cp $1  ~/Parallele/ProgComp
mpdtrace | while read ligne
do
	echo "out=$ligne";
	if [ $toto -gt $X ]; then
	echo "in=$ligne toto=$toto";
#	rm mpiuser@$ligne:Parallele/ProgComp/trexPar &
	scp trexPar mpiuser@$ligne:Parallele/ProgComp
	scp $1 mpiuser@$ligne:Parallele/ProgComp
	fi
	toto=2;
done
