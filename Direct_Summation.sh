#! /bin/bash
#
# BASH script ....

echo
echo "  Compiling the code"
gcc -Wall -g -fopenmp -DMODE_PROD Direct_Summation.c -lm -o Direct_Summation.x

if [ $? -eq 0 ]
then
  for x in $(seq 4 -1 0)
  do
    OMP_NUM_THREADS=$(echo "2^${x}" | bc) 
    echo "    Executing the code for OMP_NUM_THREADS=${OMP_NUM_THREADS}"
    # ./KE2d.x
    ./Direct_Summation.x 2>&1 | tee Direct_Summation_${OMP_NUM_THREADS}.txt
    sleep 5
  done
fi
