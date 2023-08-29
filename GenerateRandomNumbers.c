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
    // 定义生成随机数组的参数
    int array_dim1 = 2048, array_dim2 = 2048;
    float* random_array = (float*)malloc(array_dim1 * array_dim2 * sizeof(float));

    // 调用函数生成随机数组
    generate_random_array(random_array, array_dim1, array_dim2);

    // 创建并打开.dat文件
    FILE* file = fopen("random_array.dat", "wb");
    if (file == NULL) {
        perror("Failed to open file");
        return 1;
    }
    
    // 将数组内容写入文件
    for (int i = 0; i < array_dim1; i++) {
        for (int j = 0; j < array_dim2; j++) {
            fprintf(file, "%f ", random_array[i * array_dim2 + j]);
        }
        fprintf(file, "\n");
    }

    // 关闭文件
    fclose(file);
    
    // 释放数组内存
    free(random_array);

    return 0;
}
