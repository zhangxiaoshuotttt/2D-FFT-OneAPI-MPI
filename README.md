# 2D-FFT-OneAPI-MPI
This repository contains code developed for an Intel competition. It focuses on generating a 2D array of random numbers and performing a 2D FFT using MPI and the Intel OneAPI framework.

The codebase consists of five main components. To run these codes on Ubuntu 22.04 with Intel OneAPI, simply execute the `./compile_and_run.sh` script for testing purposes.

First and foremost, I want to mention that on my computer, the parallel version(`OnlyMpiFFT`) is actually slower than directly calling the serial version of the function(`FFTW3_FFT`). I suspect this is because I only utilize 12 threads and the dataset is relatively small. It is possible that the parallel version may outperform the serial version when using a larger dataset and thousands of threads.

## Generate a 2D array with dimensions of 2048x2048
The purpose of the `GenerateRandomNumbers` function is to generate a 2D array with dimensions of 2048x2048 and save it into a file. This functionality is achieved using the OneAPI library. For more detailed information, you can refer to the code in the `GenerateRandomNumbers.c` file. And here are the core codes. Within this function, I utilize the functions from the MKL library to generate a 2D array with a size of 2048x2048.
```c++
void generate_random_array(float* array, int N1, int N2)
{
    MKL_UINT seed = 0;
    VSLStreamStatePtr stream;
    vslNewStream(&stream, VSL_BRNG_MT19937, seed);
    vsRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, N1 * N2, array, 0.0, 10.0);
    vslDeleteStream(&stream);
}
```
## Using the `FFTW3_FFT` function for fast execution with OneAPI’s integrated FFTW3 library
The `FFTW3_FFT` function is a high-performance Fast Fourier Transform (FFT) implementation provided by the integrated FFTW3 library within OneAPI. By directly calling this function, you can achieve exceptionally fast execution speeds. As mentioned earlier, here are the core codes. All the necessary functions are sourced directly from the MKL library, so they can be utilized without any additional modification.
```c++
void fftw3_2D(float* Datain, fftwf_complex* Dataout, int N1, int N2)
{
    fftwf_plan stream = NULL;
    stream = fftwf_plan_dft_r2c_2d(N1, N2, Datain, Dataout, FFTW_ESTIMATE);
    fftwf_execute(stream);
    fftwf_destroy_plan(stream);
}
```
## Implementing FFT with MPI: Building Codes for Parallel Computing
Before delving into the MPI version, if you’re a newcomer like me, I highly recommend reading this [tutorial](https://mpitutorial.com/tutorials/). It was tremendously helpful to me. According to the MKL library, you can accomplish FFT using the following code, which represents the serial version. I highly recommend reading the introductory documentation for these functions as the comments are easily comprehensible. The reason I wrote these serial version codes is that, in fact, the parallel version utilizes nearly identical basic functions.

```c++
void fftMKL(float* Datain, MKL_Complex8* Dataout, int N1, int N2)
{
    MKL_LONG dims[2] = { N2, N1 };
    MKL_LONG status = 0;
    MKL_LONG rstrides[] = { 0, N1, 1 };
    MKL_LONG cstrides[] = { 0, N1 / 2 + 1, 1 };
    DFTI_DESCRIPTOR_HANDLE hand = NULL;
    char version[DFTI_VERSION_LENGTH];

    //Basic MLK functions
    DftiGetValue(0, DFTI_VERSION, version);
    status = DftiCreateDescriptor(&hand, DFTI_SINGLE, DFTI_REAL, 2, dims);
    status = DftiSetValue(hand, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    status = DftiSetValue(hand, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
    status = DftiSetValue(hand, DFTI_INPUT_STRIDES, rstrides);
    status = DftiSetValue(hand, DFTI_OUTPUT_STRIDES, cstrides);
    status = DftiCommitDescriptor(hand);
    status = DftiComputeForward(hand, Datain, Dataout);
    DftiFreeDescriptor(&hand);
    
}
```

The `OnlyMpiFFT` function is a parallel version implemented to perform 2D Fast Fourier Transform (FFT) calculations. The Intel’s instruction text provides comprehensive and detailed information about this implementation. If you are interested, I recommend reading [it](https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2023-2/fourier-transform-functions.html) to gain a better understanding of the process. If you are not familiar with MPI (Message Passing Interface), it would be beneficial to [read about it](https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2023-2/fourier-transform-functions.html) first. You can start by implementing a serial 2D FFT to become familiar with the basic functions, and then proceed to the parallel version.

**In the file `OnlyMpiFFT.c,` I have extensively commented on the parallel version of the codes, including instructions on error checking.**

Finally, the `CompareFFTResults` function is used to compare the results of the FFT calculations obtained from different versions. This can be helpful in verifying the correctness and consistency of the parallel and serial implementations of the FFT algorithm.

**These are the results obtained after running the `./compile_and_run.sh` script.**


![12核编译结果及精度对比](https://github.com/zhangxiaoshuotttt/2D-FFT-OneAPI-MPI/assets/84370046/b3e3bdeb-c8a4-4520-9c9e-07af0b948bf6)
