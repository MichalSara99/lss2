#include "lss_dense_solvers_policy.hpp"
#include "../common/lss_macros.hpp"

namespace lss_dense_solvers_policy
{

void dense_solver_qr::dense_solver_qr::solve(cusolverDnHandle_t cusolver_handle, cublasHandle_t cublas_handle,
                                             std::size_t n, const double *d_Acopy, std::size_t lda, const double *d_b,
                                             double *d_x)
{
    int bufferSize = 0;
    int bufferSize_geqrf = 0;
    int bufferSize_ormqr = 0;
    int *info = NULL;
    double *buffer = NULL;
    double *A = NULL;
    double *tau = NULL;
    int h_info = 0;
    const double one = 1.0;

    CUSOLVER_STATUS(cusolverDnDgeqrf_bufferSize(cusolver_handle, n, n, (double *)d_Acopy, lda, &bufferSize_geqrf));
    CUSOLVER_STATUS(cusolverDnDormqr_bufferSize(cusolver_handle, CUBLAS_SIDE_LEFT, CUBLAS_OP_T, n, 1, n, A, lda, NULL,
                                                d_x, n, &bufferSize_ormqr));

    bufferSize = (bufferSize_geqrf > bufferSize_ormqr) ? bufferSize_geqrf : bufferSize_ormqr;

    CUDA_ERROR(cudaMalloc(&info, sizeof(int)));
    CUDA_ERROR(cudaMalloc(&buffer, sizeof(double) * bufferSize));
    CUDA_ERROR(cudaMalloc(&A, sizeof(double) * lda * n));
    CUDA_ERROR(cudaMalloc((void **)&tau, sizeof(double) * n));

    // prepare a copy of A because getrf will overwrite A with L
    CUDA_ERROR(cudaMemcpy(A, d_Acopy, sizeof(double) * lda * n, cudaMemcpyDeviceToDevice));

    CUDA_ERROR(cudaMemset(info, 0, sizeof(int)));

    // compute QR factorization
    CUSOLVER_STATUS(cusolverDnDgeqrf(cusolver_handle, n, n, A, lda, tau, buffer, bufferSize, info));

    CUDA_ERROR(cudaMemcpy(&h_info, info, sizeof(int), cudaMemcpyDeviceToHost));

    LSS_ASSERT(h_info == 0, "LU factorization failed\n");

    CUDA_ERROR(cudaMemcpy(d_x, d_b, sizeof(double) * n, cudaMemcpyDeviceToDevice));

    // compute Q^T*b
    CUSOLVER_STATUS(cusolverDnDormqr(cusolver_handle, CUBLAS_SIDE_LEFT, CUBLAS_OP_T, n, 1, n, A, lda, tau, d_x, n,
                                     buffer, bufferSize, info));

    // x = R \ Q^T*b
    CUBLAS_STATUS(cublasDtrsm_v2(cublas_handle, CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_N,
                                 CUBLAS_DIAG_NON_UNIT, n, 1, &one, A, lda, d_x, n));
    CUDA_ERROR(cudaDeviceSynchronize());

    if (info)
    {
        CUDA_ERROR(cudaFree(info));
    }
    if (buffer)
    {
        CUDA_ERROR(cudaFree(buffer));
    }
    if (A)
    {
        CUDA_ERROR(cudaFree(A));
    }
    if (tau)
    {
        CUDA_ERROR(cudaFree(tau));
    }
}

void dense_solver_lu::dense_solver_lu::solve(cusolverDnHandle_t cusolver_handle, cublasHandle_t cublas_handle,
                                             std::size_t n, const double *d_Acopy, std::size_t lda, const double *d_b,
                                             double *d_x)
{
    int bufferSize = 0;
    int *info = NULL;
    double *buffer = NULL;
    double *A = NULL;
    int *ipiv = NULL; // pivoting sequence
    int h_info = 0;

    CUSOLVER_STATUS(cusolverDnDgetrf_bufferSize(cusolver_handle, n, n, (double *)d_Acopy, lda, &bufferSize));

    CUDA_ERROR(cudaMalloc(&info, sizeof(int)));
    CUDA_ERROR(cudaMalloc(&buffer, sizeof(double) * bufferSize));
    CUDA_ERROR(cudaMalloc(&A, sizeof(double) * lda * n));
    CUDA_ERROR(cudaMalloc(&ipiv, sizeof(int) * n));

    // prepare a copy of A because getrf will overwrite A with L
    CUDA_ERROR(cudaMemcpy(A, d_Acopy, sizeof(double) * lda * n, cudaMemcpyDeviceToDevice));
    CUDA_ERROR(cudaMemset(info, 0, sizeof(int)));

    CUSOLVER_STATUS(cusolverDnDgetrf(cusolver_handle, n, n, A, lda, buffer, ipiv, info));
    CUDA_ERROR(cudaMemcpy(&h_info, info, sizeof(int), cudaMemcpyDeviceToHost));

    LSS_ASSERT(h_info == 0, " LU factorizartion failed\n");

    CUDA_ERROR(cudaMemcpy(d_x, d_b, sizeof(double) * n, cudaMemcpyDeviceToDevice));
    CUSOLVER_STATUS(cusolverDnDgetrs(cusolver_handle, CUBLAS_OP_N, n, 1, A, lda, ipiv, d_x, n, info));
    CUDA_ERROR(cudaDeviceSynchronize());

    if (info)
    {
        CUDA_ERROR(cudaFree(info));
    }
    if (buffer)
    {
        CUDA_ERROR(cudaFree(buffer));
    }
    if (A)
    {
        CUDA_ERROR(cudaFree(A));
    }
    if (ipiv)
    {
        CUDA_ERROR(cudaFree(ipiv));
    }
}

void dense_solver_cholesky::dense_solver_cholesky::solve(cusolverDnHandle_t cusolver_handle,
                                                         cublasHandle_t cublas_handle, std::size_t n,
                                                         const double *d_Acopy, std::size_t lda, const double *d_b,
                                                         double *d_x)
{
    int bufferSize = 0;
    int *info = NULL;
    double *buffer = NULL;
    double *A = NULL;
    int h_info = 0;

    CUSOLVER_STATUS(
        cusolverDnDpotrf_bufferSize(cusolver_handle, CUBLAS_FILL_MODE_LOWER, n, (double *)d_Acopy, lda, &bufferSize));

    CUDA_ERROR(cudaMalloc(&info, sizeof(int)));
    CUDA_ERROR(cudaMalloc(&buffer, sizeof(double) * bufferSize));
    CUDA_ERROR(cudaMalloc(&A, sizeof(double) * lda * n));

    // prepare a copy of A because potrf will overwrite A with L
    CUDA_ERROR(cudaMemcpy(A, d_Acopy, sizeof(double) * lda * n, cudaMemcpyDeviceToDevice));
    CUDA_ERROR(cudaMemset(info, 0, sizeof(int)));

    CUSOLVER_STATUS(cusolverDnDpotrf(cusolver_handle, CUBLAS_FILL_MODE_LOWER, n, A, lda, buffer, bufferSize, info));

    CUDA_ERROR(cudaMemcpy(&h_info, info, sizeof(int), cudaMemcpyDeviceToHost));

    LSS_ASSERT(h_info == 0, "Cholesky factorization failed\n");

    CUDA_ERROR(cudaMemcpy(d_x, d_b, sizeof(double) * n, cudaMemcpyDeviceToDevice));

    CUSOLVER_STATUS(cusolverDnDpotrs(cusolver_handle, CUBLAS_FILL_MODE_LOWER, n, 1, A, lda, d_x, n, info));

    CUDA_ERROR(cudaDeviceSynchronize());

    if (info)
    {
        CUDA_ERROR(cudaFree(info));
    }
    if (buffer)
    {
        CUDA_ERROR(cudaFree(buffer));
    }
    if (A)
    {
        CUDA_ERROR(cudaFree(A));
    }
}

} // namespace lss_dense_solvers_policy
