#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void compare_results(const char* file1_real, const char* file1_imaginary,
                     const char* file2_real, const char* file2_imaginary)
{
    FILE* file1_real_ptr = fopen(file1_real, "rb");
    FILE* file1_imaginary_ptr = fopen(file1_imaginary, "rb");
    FILE* file2_real_ptr = fopen(file2_real, "rb");
    FILE* file2_imaginary_ptr = fopen(file2_imaginary, "rb");

    if (file1_real_ptr == NULL || file1_imaginary_ptr == NULL ||
        file2_real_ptr == NULL || file2_imaginary_ptr == NULL)
    {
        perror("Failed to open file");
        return;
    }

    int count = 0;
    double sum_real_error = 0.0;
    double sum_imaginary_error = 0.0;

    float value1_real, value1_imaginary, value2_real, value2_imaginary;
    while (fread(&value1_real, sizeof(float), 1, file1_real_ptr) == 1 &&
           fread(&value1_imaginary, sizeof(float), 1, file1_imaginary_ptr) == 1 &&
           fread(&value2_real, sizeof(float), 1, file2_real_ptr) == 1 &&
           fread(&value2_imaginary, sizeof(float), 1, file2_imaginary_ptr) == 1)
    {
        double real_error = fabs(value1_real - value2_real);
        double imaginary_error = fabs(value1_imaginary - value2_imaginary);

        sum_real_error += real_error;
        sum_imaginary_error += imaginary_error;
        count++;
    }

    if (count > 0) {
        double average_real_error = sum_real_error / count;
        double average_imaginary_error = sum_imaginary_error / count;

        printf("Overall:\n");
        printf("  Real Part: Average Error = %lf\n", average_real_error);
        printf("  Imaginary Part: Average Error = %lf\n", average_imaginary_error);

        if (average_real_error < 0.001 && average_imaginary_error < 0.001) {
            printf("The results are consistent.\n");
        }
    }

    fclose(file1_real_ptr);
    fclose(file1_imaginary_ptr);
    fclose(file2_real_ptr);
    fclose(file2_imaginary_ptr);
}

int main()
{
    const char* file1_real = "fftw3_fft_output_real_part.txt";
    const char* file1_imaginary = "fftw3_fft_output_imaginary_part.txt";
    const char* file2_real = "mkl_fft_output_real_part.txt";
    const char* file2_imaginary = "mkl_fft_output_imaginary_part.txt";

    compare_results(file1_real, file1_imaginary, file2_real, file2_imaginary);

    return 0;
}
