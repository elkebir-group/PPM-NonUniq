# VAFFP-NonUniq

Steps to generate simulations:

1. Run simulate with parameters:


```
Cover=10000

for M in {0.1,1,10}
do
    if [ ! -e simulate_output/M$M ]
    then
        mkdir simulate_output/M$M
    fi
    echo simulate_output/M$M
    for S in {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30}
    do
        ./executables/simulate -kP 10 -mut $M -k 0 -m 1 -s $S -P 1 -C $Cover > simulate_output/M$M/sim_S$S.txt
	#touch sequence_output/seq_M$M_S$S.txt
        ./executables/tree2freqs simulate_output/M${M}/sim_S${S}.txt > tree2freqs_output/freqs_M${M}_S${S}.txt
        #./executables/sequence -C $Cover simulate_output/M$M/sim_S$S.txt > sequence_output/seq_M${M}_S${S}.txt
    done
done:w

```

2.


