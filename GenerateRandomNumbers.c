#include <stdio.h>
#include <stdlib.h>
#include "mkl.h"
#include "mkl_vsl.h"
#include "mkl_service.h"

void generate_random_array(float* array, int N1, int N2)
{
    MKL_UINT seed = 0;
    VSLStreamStatePtr stream;
    vslNewStream(&stream, VSL_BRNG_MT19937, seed);
    vsRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, N1 * N2, array, 0.0, 10.0);
    vslDeleteStream(&stream);
}

int main(int argc, char *argv[])
{
    // Define the parameters for generating random arrays
    int array_dim1 = 2048, array_dim2 = 2048;
    float* random_array = (float*)malloc(array_dim1 * array_dim2 * sizeof(float));

    // Call the function to generate random arrays
    generate_random_array(random_array, array_dim1, array_dim2);

    // Create and open the .dat file
    FILE* file = fopen("random_array.dat", "wb");
    if (file == NULL) {
        perror("Failed to open file");
        return 1;
    }
    
    // Write the array content to the file
    for (int i = 0; i < array_dim1; i++) {
        for (int j = 0; j < array_dim2; j++) {
            fprintf(file, "%f ", random_array[i * array_dim2 + j]);
        }
        fprintf(file, "\n");
    }

    // Close the file
    fclose(file);
    
    // Release the memory of the array
    free(random_array);
    return 0;
}
