#!/usr/bin/env bash
# ------------------------------------------------------------------------
# SimNK version 0.1
# Project AevolStatPhys
# ------------------------------------------------------------------------
# runsAnnealed.sh / version_0.1
# ------------------------------------------------------------------------
# This script runs the executable nk_walk automatically to simulate 
# different random instances of the fitness landscape (annealed statistic).
# ------------------------------------------------------------------------
# The command line to run the script is for example:  
# Linux: ~$ bash runsAnnealed.sh
# macOS: ~$ sh runsAnnealed.sh (tested in macOS Mojave v.10.14.6)
# ------------------------------------------------------------------------
# In general: A SCRIPT TO READ THE ARGUMENTS, EXECUTE THE RUNS, SAVE THE 
# DATA, WRITE SIMULATIONS REPORTS IN README FILES AND SAVE THE SEEDS USED 
# IN THE PSEUDO RANDOM NUMBER GENERATOR. EVERYTHING CAN BE SHARED VIA 
# DROPBOX
# ------------------------------------------------------------------------
# this script samples different instances of the fitness landscape
# (the seed is random - $RANDOM - )
# ------------------------------------------------------------------------
# Nomenclature:
# ./nk_walk  -n $N -k $K -a $A -e EPI -snk SEED_NK -swlk $SEED_WLK 
#            -t $T -m M 
# N:        The genome length.
#	K:        The number of epistatic interactions (0 <= K < N).
#	A:        The alphabet size (default = 2).
#	epi:      The type of epistatic interactions ('ADJ' (default) or 'RND').
#	nk_seed:  The seed value for the landscape (default = -1).
#	wlk_seed: The seed value for the random walk (default = -1).
#	T:        The number of steps in the random walk (default = 10000).
#	M:        Mutation type (float) (0.0 <= M <= 1.0)," 
#            such that M = 0.0: point mutations only (default)," 
#                      M = 1.0: inversions only."
#            Remark: inversions of length 1 are allowed."
# ------------------------------------------------------------------------
# Created: 2019-09-01 
# Modified: 2021-02-01
# by Leonardo Trujillo (leonardo.trujillo@gmail.fr)
# ------------------------------------------------------------------------
# The Matrix has you...Follow the white rabbit. Knock, Knock, Neo

# ------------------------------------------------------------------------
# Start displaying the vintage terminal presentation
clear
echo 
echo
echo "You are welcome to our binary world!"
echo " "
sleep 0.95
clear
echo " ------------------------------------------------------"
echo " |                           /\                       |"
echo " |      SimNK       /\      /  \  /\                  |"
echo " |                 /  \_/\_/    \/  \                 |"
echo " |  0101101010100 /                  \ 1010010110011  |"           
echo " |  1010010110011/                    \0101101010100  |"
echo " |  ____________/                      \_____________ |"
echo " | Explore rugged fitness landscapes through adaptive |" 
echo " | walks to simulate evolution by natural selection.  |"
echo " |                                                    |"
echo " | version 0.1                                        |"   
echo " | 2019-2021                                          |"
echo " |                                                    |"
echo " | By Leonardo Trujillo                               |"
echo " ------------------------------------------------------"
echo
echo "> to start press ENTER and have fun doing science!!!!!!"
echo
read 
clear

# ------------------------------------------------------------------------
# Read parametes for the simulation
echo "------------------------------------------------------------"
echo " Please introduce the parameters to run the simulations."
echo " (press Ctrl C if you want to get out of our binary world.)"
echo "------------------------------------------------------------"
echo "> Genome length N:                                         |"
read N
N=$N
echo "------------------------------------------------------------"
echo "> Initial value of K assumed to be 0."
ini_K=0
echo "------------------------------------------------------------"
echo "> Maximun value of K assumed to be N."
K_max=$N-1
echo "------------------------------------------------------------"
echo "> Increments of K: assumed to be 1."
Delta_K=1 
echo "------------------------------------------------------------"
echo "> Default alphabet size A:" 
echo "2"
A=2
echo "------------------------------------------------------------"
echo "> Epistatic interactions ('ADJ' or 'RND')"
read epi
epi=$epi 
echo "------------------------------------------------------------"
echo "> shell RANDOM seed:"
#read seed
seed= #$seed

# The next part is not used in our specific experiment but can be a 
# valuable toolbox for other experiments

#echo "------------------------------------------------------------"
#echo "> Seed value for the landscape:"
#read snk
#snk=$snk
#echo "------------------------------------------------------------"
#echo "> Seed value for the random walk"
#read swlk
#swlk=$swlk 
#echo "------------------------------------------------------------"
#echo "> Iterations  (MaxIter ~ N^3):" 


#read T
T=5000 #$T 
echo "------------------------------------------------------------"
echo "> Fraction of mutations (a float number between 0.0 and 1.0)"
echo "> (0.0 for point mutations alone or 1.0 for inversion alone)"
echo ">>>In this special phase we test both ! " 
echo "------------------------------------------------------------"
echo "> Sample number of landscape instances: hard coded to 100."
#read samples
samples=100 #$samples
echo "------------------------------------------------------------"

#echo "> Cloud storage in Dropbox? "
#echo "type y if yes,  or n if not"
#share = n

echo "------------------------------------------------------------"
echo
echo "Thank you very much and press ENTER to start the simulations"
echo
echo
read
clear
echo "Here we go!"
sleep 0.9
clear 

# ------------------------------------------------------------------------
echo "---------------------------------------------------"
echo "> Simulation in progress for annealed statistics"
echo "> running $samples random hikes"
echo "> |N=$N|K=$K|A=$A|epi=$epi|snk=random|swlk=$swlk|T=$T|m='both'|samples=$samples|"
echo "> shell RANDOM seed = $seed"

now=$(date +"%T")
startTime=$now
dateLabel=$(date +"%Y_%m_%d")

echo "> Current time : $now"
echo "---------------------------------------------------"

RANDOM=$seed

# ------------------------------------------------------------------------
#mkdir Exp_"$dateLabel"
mkdir N"$N"_Annealed_"$epi" 
mkdir N"$N"_Annealed_"$epi"/DATA
#mkdir N"$N"_Annealed_"$epi"/README
mkdir N"$N"_Annealed_"$epi"/SEEDS

# ------------------------------------------------------------------------
# K evaluation cycles 
for((K=ini_K; K<=K_max; K=K+Delta_K))
do
  runLabel=N"$N"_K"$K"_T"$T"_epi"$epi"

  # ------------------------------------------------------------------------
  # Advertise the beginning of the simulation for a given K
  echo "---------------------------------------------------"
  echo "now running $runLabel"
  echo "> K = $K"
      
  # ------------------------------------------------------------------------
  # Performs the sampling
  for ((i=1;i<=samples;i++))
  do
    snk=$RANDOM
    swlk=$RANDOM
    echo "hike: $i"_"$runLabel"
    ./C++code/nk_walk -n $N -k $K  -a $A -e $epi -snk $snk  -swlk $swlk -t $T -m 1 -fullstop
    ./C++code/nk_walk -n $N -k $K  -a $A -e $epi -snk $snk  -swlk $swlk -t $T -m 0 -fullstop
  done;
done;  


    
# ------------------------------------------------------------------------
# Advertise the end of the simulation for a given K
echo "---------------------------------------------------"
echo "> Simulation finished"

now=$(date +"%T")
day=$(date +"%Y/%m/%d")

echo "> Current time : $now"
echo "> Date : $day"
echo "---------------------------------------------------"


mv final_fitness.csv final_fitness_N"$N"_"$epi".csv
 

#if [ $share == y ]
#  then 
#  tar -zcvf N"$N"Annealed.tar.gz N"$N"Annealed/
#  mv N"$N"Annealed.tar.gz ~leo/Dropbox/
#fi
 
#printf 'Working'
#for ((i = 0; i < 5; ++i)); do
#    for ((j = 0; j < 4; ++j)); do
#        printf .
#        sleep 5
#    done
#
#    printf '\b\b\b\b    \b\b\b\b'
#done
#printf '....done\n'
#$
# ------------------------------------------------------------------------
# Bye!
# ------------------------------------------------------------------------

