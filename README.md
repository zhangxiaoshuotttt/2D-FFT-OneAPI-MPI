# 2D-FFT-OneAPI-MPI
This repository contains code developed for an Intel competition. It focuses on generating a 2D array of random numbers and performing a 2D FFT using MPI and the Intel OneAPI framework.

The codebase consists of five main components. To run these codes on Ubuntu 22.04 with Intel OneAPI, simply execute the `./compile_and_run.sh` script for testing purposes.


First and foremost, I want to mention that on my computer, the parallel version(`OnlyMpiFFT`) is actually slower than directly calling the serial version of the function(`FFTW3_FFT`). I suspect this is because I only utilize 12 threads and the dataset is relatively small. It is possible that the parallel version may outperform the serial version when using a larger dataset and thousands of threads.

The purpose of the `GenerateRandomNumbers` function is to generate a 2D array with dimensions of 2048x2048 and save it into a file. This functionality is achieved using the OneAPI library. For more detailed information, you can refer to the code in the `GenerateRandomNumbers.c` file.

The `FFTW3_FFT` function is an internal integrated FFTW3 library within OneAPI. You can directly call this function, and the execution speed is still very fast.

The `OnlyMpiFFT` function is a parallel version implemented to perform 2D Fast Fourier Transform (FFT) calculations. The Intel’s instruction text provides comprehensive and detailed information about this implementation. If you are interested, I recommend reading [it](https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2023-2/fourier-transform-functions.html) to gain a better understanding of the process. If you are not familiar with MPI (Message Passing Interface), it would be beneficial to [read about it](https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2023-2/fourier-transform-functions.html) first. You can start by implementing a serial 2D FFT to become familiar with the basic functions, and then proceed to the parallel version.

Finally, the `CompareFFTResults` function is used to compare the results of the FFT calculations obtained from different versions. This can be helpful in verifying the correctness and consistency of the parallel and serial implementations of the FFT algorithm.
