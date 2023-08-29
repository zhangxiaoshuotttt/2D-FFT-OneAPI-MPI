#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "fftw3.h"
#include <time.h>
#include "mkl.h"

void fftw3_2D(float* Datain, fftwf_complex* Dataout, int N1, int N2)
{
    fftwf_plan stream = NULL;
    stream = fftwf_plan_dft_r2c_2d(N1, N2, Datain, Dataout, FFTW_ESTIMATE);
    fftwf_execute(stream);
    fftwf_destroy_plan(stream);
}

float* My_read_data_file_2d(const char* filename, int m, int n)
{
    // open file
    FILE* file = fopen(filename, "r");
    if (file == NULL) {
        perror("Failed to open file");
        return NULL;
    }

    // allocate memory for array
    float* array = (float*)malloc(m * n * sizeof(float));
    if (array == NULL) {
        perror("Failed to allocate memory");
        fclose(file);
        return NULL;
    }

    // read data from file to array
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            if (fscanf(file, "%f", &(array[i * n + j])) != 1) {
                perror("Failed to read data from file");
                fclose(file);
                free(array);
                return NULL;
            }
        }
    }

    // close file
    fclose(file);

    return array;
}

void My_save_data_file_2d(const char* filename, fftwf_complex* array, long size)
{
    char real_part_file[256];
    char imaginary_part_file[256];
    sprintf(real_part_file, "%s_real_part.txt", filename);
    sprintf(imaginary_part_file, "%s_imaginary_part.txt", filename);

    FILE* real_part_file_ptr = fopen(real_part_file, "w");
    FILE* imaginary_part_file_ptr = fopen(imaginary_part_file, "w");

    if (real_part_file_ptr == NULL || imaginary_part_file_ptr == NULL) {
        perror("Failed to open file");
        return;
    }

    // 将实部和虚部分别写入到两个不同的文件中
    for (long i = 0; i < size; i++) {
        float real_part = fabsf(array[i][0]);
        float imaginary_part = fabsf(array[i][1]);

        fprintf(real_part_file_ptr, "%f\n", real_part);
        fprintf(imaginary_part_file_ptr, "%f\n", imaginary_part);
    }

    fclose(real_part_file_ptr);
    fclose(imaginary_part_file_ptr);
}

int main(int argc, char *argv[])
{
    int array_dim1 = 2048, array_dim2 = 2048;
    float* random_array = My_read_data_file_2d(argv[1], array_dim1, array_dim2);
    if (random_array == NULL) {
        return 1;
    }
    
    float* fftw3_input = (float*)fftwf_malloc(array_dim1 * array_dim2 * sizeof(float));
    fftwf_complex* fftw3_output = (fftwf_complex*)fftwf_malloc((array_dim1)*array_dim2 * sizeof(fftwf_complex));
    for (int j = 0; j < array_dim1 * array_dim2; j++)
    {
        fftw3_input[j] = random_array[j];
    }
    

    // perform FFT
    clock_t start_time_FFTW = clock();

    fftw3_2D(fftw3_input, fftw3_output, array_dim1, array_dim2);

    clock_t end_time_FFTW = clock();
    float total_time_FFTW = (float)(end_time_FFTW - start_time_FFTW) / CLOCKS_PER_SEC;
    printf("PerformFFT_fftw total time: %f s \n", total_time_FFTW);
    
    My_save_data_file_2d("fftw3_fft_output",fftw3_output,array_dim1*array_dim2);


    // free memory
    free(random_array);
    fftwf_free(fftw3_output);

    return 0;
}
