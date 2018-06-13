#!/bin/bash
if [ ! -e simulate_output ]
then
    mkdir simulate_output
fi

for M in {0.1,1,10}
do
    if [ ! -e simulate_output/M$M ]
    then
        mkdir simulate_output/M$M
    fi
    #./executables/simulate -kP 1 -mut $M -k 0 -m 1 -N 20 -P 1 -p 0 -C 1000000 -o simulate_output/M$M/
done

if [ ! -e mix_output ]
then
    mkdir mix_output
fi

for M in {0.1,1,10}
do
    for f in simulate_output/M$M/*.tree
    do
        for k in {1,2,5,10}
        do
            seed=$(basename $f .tree | sed -e "s/T_seed\([0-9]*\)/\1/g")
            ff=mix_output/M${M}_S${seed}_k${k}.tree
            ./executables/mix -p -k $k -s 0 $f > $ff
        done
    done
done

if [ ! -e tree2freqs_output ]
then
    mkdir tree2freqs_output
fi

for f in mix_output/*.tree
do
    ./executables/tree2freqs $f > tree2freqs_output/$(basename $f .tree).tsv
done
