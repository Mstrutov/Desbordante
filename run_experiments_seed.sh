#!/bin/bash

declare -a seeds=(1000 2000 10000)
dataset=CIPublicHighway

cd build/target

declare experiment_name=01_16_seed
mkdir $experiment_name
for j in "${seeds[@]}"
do
	#echo $j
	#cat 01_07/$j/results_${datasets[i]}.txt | grep 'TIME' > 01_07/$j/results_${datasets[i]}.txt
	touch $experiment_name/results_$dataset.txt
	#echo --algo='pyro' --data="${datasets[i]}".csv --sep="${separators[i]}" --hasHeader="${header_presence[i]}" --seed=$j 
	./fdtester_run --algo='pyro' --data=$dataset.csv --sep=',' --seed=$j >> $experiment_name/results_$dataset.txt
done

declare trunc_name=trunc_$experiment_name
cat $experiment_name/results_$dataset.txt | grep 'TIME' > $experiment_name/results_trunc_$dataset.txt
