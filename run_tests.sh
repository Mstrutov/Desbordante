#!/bin/bash
declare -a datasets=('adult' 'breast_cancer' 'CIPublicHighway' 'EpicMeds' 'EpicVitals' 'iowa1kk' 'LegacyPayors' 'neighbors100k')
declare -a separators=(';' ',' ',' '|' '|' ',' '|' ',')
declare -a header_presence=('false' 'true' 'true' 'true' 'true' 'true' 'true' 'true')

declare -a seeds=(10000)

cd build/target

declare experiment_name=02_08_nocaching
mkdir $experiment_name
for j in "${seeds[@]}"
do
mkdir $experiment_name/$j
for i in $(seq 0 7)
do
	#cat 01_07/$j/results_${datasets[i]}.txt | grep 'TIME' > 01_07/$j/results_${datasets[i]}.txt
	for k in $(seq 1 1)
	do
		touch $experiment_name/$j/results_${datasets[i]}.txt
		./Desbordante_run --task=fd --algo='pyro' --data="${datasets[i]}".csv --sep="${separators[i]}" --has_header="${header_presence[i]}" --seed=$j | grep 'TIME' >> $experiment_name/$j/results_${datasets[i]}.txt
	done
done
done

# перепрогнать на EpicVitals
