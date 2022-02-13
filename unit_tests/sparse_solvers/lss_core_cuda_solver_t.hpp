#if !defined(_LSS_CORE_CUDA_SOLVE_T_HPP_)
#define _LSS_CORE_CUDA_SOLVE_T_HPP_

// ON DEVICE:
extern void device_sparse_qr_test();

extern void device_sparse_qr_pointer_test();

extern void test_device_sparse_qr();

// ON HOST:
extern void host_sparse_qr_test();

extern void host_sparse_qr_pointer_test();

extern void test_host_sparse_qr_test();

// ON HOST:
extern void host_bvp_dirichlet_bc_qr_test();

extern void host_bvp_dirichlet_bc_lu_test();

extern void host_bvp_dirichlet_bc_cholesky_test();

extern void test_dirichlet_bc_bvp_on_host();

// ON DEVICE:
extern void device_bvp_dirichlet_bc_qr_test();

extern void device_bvp_dirichlet_bc_cholesky_test();

extern void test_dirichlet_bc_bvp_on_device();

// ON HOST:
extern void host_bvp_robin_bc_qr_test();

extern void host_bvp_robin_bc_lu_test();

extern void host_bvp_robin_bc_cholesky_test();

extern void test_robin_bc_bvp_on_host();

// ON DEVICE:
extern void device_bvp_robin_bc_qr_test();

extern void device_bvp_robin_bc_cholesky_test();

extern void test_robin_bc_bvp_on_device();

#endif ///_LSS_CORE_CUDA_SOLVE_T_HPP_
