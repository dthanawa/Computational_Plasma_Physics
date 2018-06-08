all: PIC_mkl.exe

PIC_mkl.exe:
 icc -I${MKLROOT}/include -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -DMODE_DEBUG PIC_mkl.c -o PIC_mkl.exe -std=c99 -lm
