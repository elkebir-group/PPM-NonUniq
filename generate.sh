#!/bin/bash

#Cover=10000
#
#for M in {0.1,1,10}
#do
#    if [ ! -e simulate_output/M$M ]
#    then
#        mkdir simulate_output/M$M
#    fi
#    echo simulate_output/M$M
#    for S in {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30}
#    do
#        ./executables/simulate -kP 10 -mut $M -k 0 -m 1 -s $S -P 1 -C $Cover > simulate_output/M$M/sim_S$S.txt
#	#touch sequence_output/seq_M$M_S$S.txt
#        ./executables/tree2freqs simulate_output/M${M}/sim_S${S}.txt > tree2freqs_output/freqs_M${M}_S${S}.txt
#        #./executables/sequence -C $Cover simulate_output/M$M/sim_S$S.txt > sequence_output/seq_M${M}_S${S}.txt
#    done
#done

# run R code to convert freqs to reads

for C in {-1,50,100,1000}
do
    for P in 1
    do
        for k in {2,5,10}
        do
            for M in {0.1,1,10}
            do
                for S in {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30}
                do
                    if [ -s downsample_input/freqs_M${M}_S${S}.txt ]
                    then
                        ./executables/downsample -P $P -k $k -C $C downsample_input/freqs_M${M}_S${S}.txt > downsample_output/downsampled_C${C}_P${P}_k${k}_M${M}_S${S}.txt
                    fi
                    if [ -s downsample_output/downsampled_C${C}_P${P}_k${k}_M${M}_S${S}.txt ]
                    then
                        ./executables/reads2freqs downsample_output/downsampled_C${C}_P${P}_k${k}_M${M}_S${S}.txt > reads2freqs_output/freq_C${C}_P${P}_k${k}_M${M}_S${S}.txt
                    fi
                done
            done
        done
    done
done

# run r code to convert downsampled reads to freqs for spruce input

#for C in {-1,50,100,1000}
#do
#    for P in 1
#    do
#        for k in {2,5,10}
#        do
#            for M in {0.1,1,10}
#            do
#                for S in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30}
#                do
#                    #./executables/enumerate spruce_input/downsampled_C${C}_P${P}_k${k}_M${M}_S${S}.txt > enumerate_output/sim_C${C}_P${P}_k${k}_M${M}_S${S}.txt
#                    #./executables/parse enumerate_output/sim_C${C}_P${P}_k${k}_M${M}_S${S}.txt
#                done
#            done
#        done
#    done
#done 

#./parse ../result/n_15_noisy/sims_r0_m10_n15_c1000.res 

#./recall ../result/n_15_noisy/sims_r0_m10_n15_c1000.res ../data/sims/n_15_noisy/sims_r0_m10_n15_c1000.true >> trees/sims_r0_m10_n15_c1000_recall.csv
# can spruce use oncosim tree as input for measuring recall?

