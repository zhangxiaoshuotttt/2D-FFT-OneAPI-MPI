#!/bin/bash

# Compile GenerateRandomNumbers.c
icx -o GenerateRandomNumbers GenerateRandomNumbers.c -lmkl_rt

# Run GenerateRandomNumbers to generate random arrays
./GenerateRandomNumbers

# Compile FFTW3_FFT.c
icx -o FFTW3_FFT FFTW3_FFT.c -lmkl_rt

# Perform FFTW3_FFT on random_array.dat for 10 runs and record the execution time
times=10

FFTWtotal_runtime=0
for ((i=1; i<=times; i++))
do
  output=$(./FFTW3_FFT random_array.dat)
  runtime=$(echo $output | grep -oE 'PerformFFT_fftw total time: ([0-9]+\.[0-9]+)' | grep -oE '[0-9]+\.[0-9]+')

  if [ -n "$runtime" ]; then
    FFTWtotal_runtime=$(awk "BEGIN{print $FFTWtotal_runtime + $runtime}")
    echo "FFTW3_FFT running time for the $i-th run is $runtime seconds"
  else
    echo "Failed to extract the runtime for the $i-th run of FFTW3_FFT"
  fi
done
echo "The total runtime of all FFTW3_FFT runs is $FFTWtotal_runtime seconds"

# Compile OnlyMpiFFT.c
mpiicc OnlyMpiFFT.c -o OnlyMpiFFT -L${MKLROOT}/lib/intel64 -lmkl_cdft_core -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -lmkl_blacs_intelmpi_ilp64 -liomp5 -lpthread -lm -ldl

# Execute 'mpirun -np 12 OnlyMpiFFT random_array.dat -DMKL_ILP64 -I”${MKLROOT}/include"’ 10 times and record the execution time.
MPIruntimes=()
for ((i=1; i<=times; i++))
do
  MPIoutput=$(mpirun -np 12 OnlyMpiFFT random_array.dat -DMKL_ILP64 -I"${MKLROOT}/include")
  runtime=$(echo "$MPIoutput" | grep -oE 'All processor for FFT run ([0-9]+\.[0-9]+)' | grep -oE '[0-9]+\.[0-9]+')

  if [ -n "$runtime" ]; then
    MPIruntimes+=($runtime)
    echo "MPI runtime for the $i-th run is $runtime seconds"
  else
    echo "Unable to retrieve the runtime for the $i-th MPI run"
  fi
done

OnlyMpiFFT_total_time=0
for runtime in "${MPIruntimes[@]}"
do
    OnlyMpiFFT_total_time=$(awk "BEGIN{print $OnlyMpiFFT_total_time + $runtime}")
done

echo "The total runtime of all MPI runs is $OnlyMpiFFT_total_time seconds"

# Compile CompareFFTResults.c
icx -o CompareFFTResults CompareFFTResults.c -lmkl_rt

# Perform CompareFFTResults
./CompareFFTResults
