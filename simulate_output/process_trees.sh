#!/bin/bash
#for M in {0.1,0.2,0.4}
#do
#    for f in M$M/*.tree
#    do
#        S=$(echo $(basename $f .tree) | sed -e s/T_seed//g)
#        python process_tree.py $f > M${M}_S${S}.tree
#    done
#done

#for n in {3,5,7,9,11,13}
for n in {15,17,19}
do
    for f in n$n/*.tree
    do
        S=$(echo $(basename $f .tree) | sed -e s/T_seed//g)
        python process_tree.py $f > n${n}_S${S}.tree
    done
done
