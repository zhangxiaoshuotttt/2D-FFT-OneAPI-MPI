#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "mkl_service.h"
#include "mkl_dfti.h"
#include "mkl.h"
#include "mkl_vsl.h"
#include "fftw3.h"
#include <time.h>
#include "mpi.h"
#include "mkl_cdft.h"
#define PREC DFTI_SINGLE
#define ADVANCED_DATA_PRINT 1
#define ACCURACY_PRINT 1
#define LEGEND_PRINT 1

#define SINGLE_EPS 1.0E-6

typedef struct {
    float re;
    float im;
} mkl_float_complex;
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
void My_save_data_file_2d(const char* filename, MKL_Complex8* array, int m, int n)
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
    for (long i = 0; i < m * n; i++) {
        float real_part = fabsf(array[i].real);
        float imaginary_part = fabsf(array[i].imag);

        fprintf(real_part_file_ptr, "%f\n", real_part);
        fprintf(imaginary_part_file_ptr, "%f\n", imaginary_part);
    }

    fclose(real_part_file_ptr);
    fclose(imaginary_part_file_ptr);
}


int MKL_Data(MPI_Comm Comm,int RootRank,int ElementSize,long Dim,
                  MKL_LONG Lengths[],void *global,long nx,long start_x,void *local,
                  int Flag)
{

    int nProc,nRank,MPI_err,i,fd;
    int tmp[2],*counts,*displs,*buf;
    MPI_Request req;
    MPI_Status status;
    counts = displs = buf = NULL;

    MPI_err=MPI_Comm_rank(Comm,&nRank);
    if (MPI_err != MPI_SUCCESS) return MPI_err;

    if (nRank == RootRank) {
        MPI_err=MPI_Comm_size(Comm,&nProc);
        if (MPI_err != MPI_SUCCESS) return MPI_err;
        counts=(int*)mkl_malloc(nProc*sizeof(int), 64);
        displs=(int*)mkl_malloc(nProc*sizeof(int), 64);
        buf=(int*)mkl_malloc(2*nProc*sizeof(int), 64);
        if (counts == NULL || displs == NULL || buf == NULL) {
            printf("Not enough memory\n");
            return 1;
        }
    }

    fd=ElementSize;
    for (i=1;i<(int)Dim;i++) fd*=(int)Lengths[i];

    tmp[0]=(int)nx*fd;
    tmp[1]=(int)start_x*fd;

    MPI_err=MPI_Gather(tmp,2,MPI_INT,buf,2,MPI_INT,RootRank,Comm);
    if (MPI_err != MPI_SUCCESS) return MPI_err;

    if (nRank == RootRank) {
        for (i=0;i<nProc;i++) {
            counts[i]=buf[2*i];
            displs[i]=buf[2*i+1];
        }
    }

    if (Flag == 0) {
        MPI_err = MPI_Irecv(local,tmp[0],MPI_BYTE, RootRank, 123,Comm,&req);
        if (MPI_err != MPI_SUCCESS) return MPI_err;

        if (nRank == RootRank) {
            for (i=0;i<nProc;i++) {
                MPI_err = MPI_Send(((char*)global)+displs[i],counts[i],MPI_BYTE,i,123,Comm);
                if (MPI_err != MPI_SUCCESS) return MPI_err;
            }
        }
        MPI_Wait(&req,&status);
        if (MPI_err != MPI_SUCCESS) return MPI_err;
    }
    else {
        MPI_err = MPI_Gatherv(local,tmp[0],MPI_BYTE,global,counts,displs,MPI_BYTE,RootRank,Comm);
        if (MPI_err != MPI_SUCCESS) return MPI_err;
    }

    if (nRank == RootRank) {
        mkl_free(buf);
        mkl_free(displs);
        mkl_free(counts);
    }
    return 0;
}

int MKL_ScatterData(MPI_Comm Comm,int RootRank,int ElementSize,long Dim,
                         MKL_LONG Lengths[],void *global_in,long nx,long start_x,
                         void *local_in)
{
    return MKL_Data(Comm,RootRank,ElementSize,Dim,Lengths,
                         global_in,nx,start_x,local_in,0);
}

int MKL_GatherData(MPI_Comm Comm,int RootRank,int ElementSize,long Dim,
                        MKL_LONG Lengths[],void *global_out,long nx,long start_x,
                        void *local_out)
{
    return MKL_Data(Comm,RootRank,ElementSize,Dim,Lengths,
                         global_out,nx,start_x,local_out,1);
}
void status_print(int rank, long status)
{
    printf("Rank %d: Error message: %s\n", rank, DftiErrorMessage(status));
    printf("Rank %d: Error status = %d\n", rank, (int)status);
}

int main(int argc, char *argv[])
{

    int array_dim1 = 2048, array_dim2 = 2048;
    int out_dim1 = 1024, out_dim2 = 1024;
    float* random_array = My_read_data_file_2d(argv[1], array_dim1, array_dim2);
    if (random_array == NULL) {
        return 1;
    }
    
    //We need initial parameters
    float *mkl_input_real;
    MKL_Complex8* mkl_output_cmplx;
    int local_failure = 0;
    int global_failure = 0;
    mkl_input_real = NULL;
    mkl_output_cmplx=NULL;
    MKL_LONG Status = DFTI_NO_ERROR;
    MKL_LONG lengths[2];
    
    //parameters for parralell
    DFTI_DESCRIPTOR_DM_HANDLE desc;
    long RootRank,ElementSize;
    MKL_LONG nx,nx_out,start_x,start_x_out,size;
    
    //parameters for MPI
    int MPI_err,err;
    int MPI_nProc;
    int MPI_Rank;

    /**********************************************************************************
     1. Initiate MPI 
    **********************************************************************************/
    MPI_err = MPI_Init(&argc, &argv);

    if (MPI_err != MPI_SUCCESS) {
        printf(" MPI initialization error\n");
        printf(" CHEKE FAILED\n");
        return 1;
    }
    
    MPI_Comm_size(MPI_COMM_WORLD, &MPI_nProc);
    MPI_Comm_rank(MPI_COMM_WORLD, &MPI_Rank);
    /*
    if (MPI_Rank == 0)
        printf( " Program is running on %d processes\n", MPI_nProc);
    */
    if (MPI_Rank == 0) err = 0;
    MPI_Bcast(&err,1,MPI_INT,0,MPI_COMM_WORLD);
    if (err != 0) {
        global_failure++;
        goto CLOSE_MPI;
    }
    
    MPI_Bcast(&out_dim1,1,MPI_LONG_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&out_dim2,1,MPI_LONG_INT,0,MPI_COMM_WORLD);

    lengths[0] = out_dim1;
    lengths[1] = out_dim2;
    

    /*********************************************************************************
     2. Allocate memory for input real numbers and output complex numbers
    **********************************************************************************/

    //For mkl based on mpi, we need define twice for total and partial
    float* total_in = (float*)mkl_malloc(array_dim2 * array_dim1 * sizeof(float), 64);
    for (int j = 0; j < (array_dim1) * array_dim2; j++)
    {
        total_in[j] = random_array[j];
    }
    
    MKL_Complex8* total_exp = (MKL_Complex8*)mkl_malloc(out_dim1 * out_dim2 * sizeof(MKL_Complex8), 64);
    if (total_in == NULL || total_exp == NULL) {
        printf(" Rank %d: Not enough memory\n", MPI_Rank);
        local_failure++;
    }
    MPI_Allreduce(&local_failure, &global_failure, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (global_failure != 0) goto FREE_MEM;
    
    /*********************************************************************************
     3. Allocate memory for the descriptor
    *********************************************************************************/
    
    Status = DftiCreateDescriptorDM(MPI_COMM_WORLD,&desc,PREC,DFTI_REAL,2,lengths);
    /*
    if (ADVANCED_DATA_PRINT && (MPI_Rank == 0))
        printf("Create=%ld\n",(long)Status);
    if (!DftiErrorClass(Status, DFTI_NO_ERROR)) local_failure++;
    MPI_Allreduce(&local_failure, &global_failure, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (global_failure != 0) goto FREE_DESCRIPTOR;
    
    //printf("We have allocated memory for the descriotor \n");
    */
    /*******************************************************************************
     4. Get some values of configuration parameters to make sure paralell programe
        can run successfully
    ********************************************************************************/
    Status = DftiGetValueDM(desc,CDFT_LOCAL_SIZE,&size);
    /*
    if (ADVANCED_DATA_PRINT && (MPI_Rank == 0))
        printf("Get=%ld,size=%ld\n",(long)Status,(long)size);
    if (!DftiErrorClass(Status, DFTI_NO_ERROR)) local_failure++;
    MPI_Allreduce(&local_failure, &global_failure, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (global_failure != 0) goto FREE_DESCRIPTOR;
    */
    Status = DftiGetValueDM(desc,CDFT_LOCAL_NX,&nx);
    /*
    if (ADVANCED_DATA_PRINT && (MPI_Rank == 0))
        printf("Get=%ld,nx=%ld\n",(long)Status,(long)nx);
    if (!DftiErrorClass(Status, DFTI_NO_ERROR)) local_failure++;
    MPI_Allreduce(&local_failure, &global_failure, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (global_failure != 0) goto FREE_DESCRIPTOR;
    */
    Status = DftiGetValueDM(desc,CDFT_LOCAL_X_START,&start_x);
    /*
    if (ADVANCED_DATA_PRINT && (MPI_Rank == 0))
        printf("Get=%ld,start_x=%ld\n",(long)Status,(long)start_x);
    if (!DftiErrorClass(Status, DFTI_NO_ERROR)) local_failure++;
    MPI_Allreduce(&local_failure, &global_failure, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (global_failure != 0) goto FREE_DESCRIPTOR;
    */
    Status = DftiGetValueDM(desc,CDFT_LOCAL_OUT_NX,&nx_out);
    /*
    if (ADVANCED_DATA_PRINT && (MPI_Rank == 0))
        printf("Get=%ld,nx_out=%ld\n",(long)Status,(long)nx_out);
    if (!DftiErrorClass(Status, DFTI_NO_ERROR)) local_failure++;
    MPI_Allreduce(&local_failure, &global_failure, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (global_failure != 0) goto FREE_DESCRIPTOR;
    */
    Status = DftiGetValueDM(desc,CDFT_LOCAL_OUT_X_START,&start_x_out);
    /*
    if (ADVANCED_DATA_PRINT && (MPI_Rank == 0))
        printf("Get=%ld,start_x_out=%ld\n",(long)Status,(long)start_x_out);
    if (!DftiErrorClass(Status, DFTI_NO_ERROR)) local_failure++;
    MPI_Allreduce(&local_failure, &global_failure, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (global_failure != 0) goto FREE_DESCRIPTOR;
    */
    /*Here we define the local parameter for paralell programes*/
    mkl_input_real=(float*)mkl_malloc(size*sizeof(float), 64);
    mkl_output_cmplx=(MKL_Complex8*)mkl_malloc(size*sizeof(mkl_float_complex), 64);
    if (mkl_input_real == NULL || mkl_output_cmplx == NULL) {
        printf("Rank %d: Not enough memory\n", MPI_Rank);
        local_failure++;
    }
    MPI_Allreduce(&local_failure, &global_failure, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (global_failure != 0) goto FREE_DESCRIPTOR;
 
 
 
    /***********************************************************************************
     5. Set values of configuration parameters
    ***********************************************************************************/
    Status = DftiSetValueDM(desc,CDFT_WORKSPACE,mkl_output_cmplx);
    /*
    if (ADVANCED_DATA_PRINT && (MPI_Rank == 0))
        printf("Set=%ld,pointer=%p\n",(long)Status,mkl_output_cmplx);
    if (!DftiErrorClass(Status, DFTI_NO_ERROR)) local_failure++;
    MPI_Allreduce(&local_failure, &global_failure, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (global_failure != 0) goto FREE_DESCRIPTOR;
    */

    /***************************************************************************************************
     6. Perform all initialization for the actual FFT computation
    ***************************************************************************************************/
    Status = DftiCommitDescriptorDM(desc);
    /*
    if (ADVANCED_DATA_PRINT && (MPI_Rank == 0))
        printf("Commit=%ld\n",(long)Status);
    if (!DftiErrorClass(Status, DFTI_NO_ERROR)) local_failure++;
    MPI_Allreduce(&local_failure, &global_failure, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (global_failure != 0) goto FREE_DESCRIPTOR;
    */
    /***************************************************************************************************
     7. Create arrays for local parts of input and output data
    ***************************************************************************************************/
    RootRank=0;
    ElementSize=sizeof(float);
    float start_time = MPI_Wtime();
    Status = MKL_ScatterData(MPI_COMM_WORLD,RootRank,ElementSize,2,lengths,total_in,nx,start_x,mkl_input_real);   
    /*
    if (ADVANCED_DATA_PRINT && (MPI_Rank == 0))
        printf("Scatter=%ld\n",(long)Status);
    if (!DftiErrorClass(Status, DFTI_NO_ERROR)) local_failure++;
    MPI_Allreduce(&local_failure, &global_failure, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (global_failure != 0) goto FREE_DESCRIPTOR;
    */
    /**************************************************************************************************
     8. Compute the transform by calling DftiComputeForwardDM
    ***************************************************************************************************/
    
    Status = DftiComputeForwardDM(desc,mkl_input_real,mkl_output_cmplx);
    /*
    if (ADVANCED_DATA_PRINT && (MPI_Rank == 0))
        printf("ComputeForward=%ld\n",(long)Status);
    if (!DftiErrorClass(Status, DFTI_NO_ERROR)) local_failure++;
    MPI_Allreduce(&local_failure, &global_failure, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (global_failure != 0) goto FREE_DESCRIPTOR;
    */
    /************************************************************************************************
     9. Gather data among processors
    *************************************************************************************************/
    Status = MKL_GatherData(MPI_COMM_WORLD,RootRank,ElementSize,2,lengths,total_exp,nx,start_x,mkl_output_cmplx);    
    /*
    if (ADVANCED_DATA_PRINT && (MPI_Rank == 0))
        printf("Gather=%ld\n",(long)Status);
    if (!DftiErrorClass(Status, DFTI_NO_ERROR)) local_failure++;
    MPI_Allreduce(&local_failure, &global_failure, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (global_failure != 0) goto FREE_DESCRIPTOR;
     */
     /****************************output data*******************************************************/
    float end_time = MPI_Wtime();
    float runtime = end_time - start_time;
    float recv_runtime, total_runtime, average_time;
    
    MPI_Reduce(&runtime,&recv_runtime,1,MPI_FLOAT,MPI_SUM,0,MPI_COMM_WORLD);
    if(MPI_Rank==0){
        total_runtime =recv_runtime;
        average_time =total_runtime/12 ;
        printf("All processor for FFT run %f s \n",average_time);
        }
    
    My_save_data_file_2d("mkl_fft_output", total_exp, out_dim1, out_dim2);

FREE_DESCRIPTOR:

    /*
    **  Check status of DFTI functions
    */
    if (!DftiErrorClass(Status, DFTI_NO_ERROR)) {
        status_print(MPI_Rank, Status);
    }

    /*************************************************************************************************
     10. Release memory allocated for a descriptor)
    *************************************************************************************************/
    Status = DftiFreeDescriptorDM(&desc);
    /*
    if (ADVANCED_DATA_PRINT && (MPI_Rank == 0))
        printf("FreeDescriptor=%ld\n",(long)Status);
    if (! DftiErrorClass(Status, DFTI_NO_ERROR)) {
        status_print(MPI_Rank, Status);
        local_failure++;
    }
    MPI_Allreduce(&local_failure, &global_failure, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    */
FREE_MEM:

    /********************************************
       Deallocate memory for dynamic arrays
    *******************************************/

    if (total_in  != NULL) mkl_free(total_in);
    if (total_exp != NULL) mkl_free(total_exp);
    if (mkl_input_real != NULL) mkl_free(mkl_input_real);
    if (mkl_output_cmplx  != NULL) mkl_free(mkl_output_cmplx);
    //printf("We have freed the memory");

CLOSE_MPI:

    /*************************************************************************************************
     11. Finalize MPI
    *************************************************************************************************/
    MPI_Finalize();

    if (global_failure != 0) {
        if (MPI_Rank == 0) printf(" TEST FAILED\n");
        return 1;
    }

    // if (MPI_Rank == 0) printf(" END OF FFT\n");

}
