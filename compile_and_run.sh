#!/bin/bash

# 编译 GenerateRandomNumbers.c
icx -o GenerateRandomNumbers GenerateRandomNumbers.c -lmkl_rt

# 运行 GenerateRandomNumbers 以生成随机数组
./GenerateRandomNumbers

# 编译 FFTW3_FFT.c
icx -o FFTW3_FFT FFTW3_FFT.c -lmkl_rt

# 进行 FFTW3_FFT random_array.dat 运行10次，记录运行时间
times=10

FFTWtotal_runtime=0
for ((i=1; i<=times; i++))
do
  output=$(./FFTW3_FFT random_array.dat)
  runtime=$(echo $output | grep -oE 'PerformFFT_fftw total time: ([0-9]+\.[0-9]+)' | grep -oE '[0-9]+\.[0-9]+')

  if [ -n "$runtime" ]; then
    FFTWtotal_runtime=$(awk "BEGIN{print $FFTWtotal_runtime + $runtime}")
    echo "FFTW3_FFT 第 $i 次运行时间为 $runtime 秒"
  else
    echo "FFTW3_FFT 第 $i 次未能提取到运行时间"
  fi
done

echo "FFTW3_FFT 所有运行时间的总和为 $FFTWtotal_runtime 秒"

# 编译 OnlyMpiFFT.c
mpiicc OnlyMpiFFT.c -o OnlyMpiFFT -L${MKLROOT}/lib/intel64 -lmkl_cdft_core -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -lmkl_blacs_intelmpi_ilp64 -liomp5 -lpthread -lm -ldl

# 进行 mpirun -np 12 OnlyMpiFFT random_array.dat -DMKL_ILP64 -I"${MKLROOT}/include" 运行10次，记录运行时间
MPIruntimes=()
for ((i=1; i<=times; i++))
do
  MPIoutput=$(mpirun -np 12 OnlyMpiFFT random_array.dat -DMKL_ILP64 -I"${MKLROOT}/include")
  runtime=$(echo "$MPIoutput" | grep -oE 'All processor for FFT run ([0-9]+\.[0-9]+)' | grep -oE '[0-9]+\.[0-9]+')

  if [ -n "$runtime" ]; then
    MPIruntimes+=($runtime)
    echo "MPI 第 $i 次运行时间为 $runtime 秒"
  else
    echo "MPI 第 $i 次未能提取到运行时间"
  fi
done

OnlyMpiFFT_total_time=0
for runtime in "${MPIruntimes[@]}"
do
    OnlyMpiFFT_total_time=$(awk "BEGIN{print $OnlyMpiFFT_total_time + $runtime}")
done

echo "MPI 所有运行时间的总和为 $OnlyMpiFFT_total_time 秒"

# 编译 CompareFFTResults.c
icx -o CompareFFTResults CompareFFTResults.c -lmkl_rt

# 运行 CompareFFTResults
./CompareFFTResults
