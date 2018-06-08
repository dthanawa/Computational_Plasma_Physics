#! /bin/bash
#
# BASH script ....

echo
echo "  Compiling the code"
icc -I${MKLROOT}/include -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -DMODE_DEBUG PIC_mkl.c -o PIC_mkl.exe -std=c99 -lm

if [ $? -eq 0 ]
then
  ./PIC_mkl.exe 2>&1 | tee PIC_mkl.txt
fi

# Delete core dumps
# Core dumps have the file name format of core.#####
# Refer to Tip #14 in UN5390 (Spring 2018)
echo "  Deleting core dumps"
for x in $(find . -print | grep -E "core.[0-9]+$")
do
  echo "    ${x}"
  rm -f ${x}
done
