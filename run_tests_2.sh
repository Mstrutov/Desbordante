#!/bin/bash
declare -a datasets=('adult' 'breast_cancer' 'CIPublicHighway' 'EpicMeds' 'EpicVitals' 'iowa1kk' 'LegacyPayors' 'neighbors100k')
declare -a separators=(';' ',' ',' '|' '|' ',' '|' ',')
declare -a header_presence=('false' 'true' 'true' 'true' 'true' 'true' 'true' 'true')

declare -a seeds=(1000 2000 10000)

cd build/target

declare experiment_name=01_08
declare dest=trunc_$experiment_name
mkdir $dest
for j in "${seeds[@]}"
do
mkdir $dest/$j
for i in $(seq 0 7)
do
	cat $experiment_name/$j/results_${datasets[i]}.txt | grep 'TIME' > $dest/$j/results_${datasets[i]}.txt
done
done

# перепрогнать на EpicVitals
