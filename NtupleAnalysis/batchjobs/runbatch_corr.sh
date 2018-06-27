#!/bin/bash

if [[ $# -eq 0 ]] ; then
    echo 'please give 13d, 13e, or 13f as an argument'
    exit 0
fi

date

#FIXME: call Ivan's script here
module add ROOT/6.08.00

#for p in {0,4,6}
for p in 4
do
    echo "calling paired Injector with file $1.root and track Gev $p"
    ./../pair_gale_shapley/paired_injector $1.root $p
    for i in {0..280..20} #Mix 300 events
    do
	mix_min=$i
	mix_max="$((i + 19))"
	#sbatch -p shared-chos -t 16:00:00 runCorr.sh $1 $mix_min $mix_max $p
	./runCorr.sh $1 $mix_min $mix_max $p
    echo "$mix_min $mix_max $1"
    done
done