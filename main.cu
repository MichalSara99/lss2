#include<iostream>
#include<string>

// TESTS FOR UTILITY:
#include"unit_tests/common/lss_numerical_integration_t.hpp"
#include"unit_tests/common/lss_heston_price_t.hpp"
#include"unit_tests/common/lss_xml_t.hpp"
// TESTS FOR CONTAINERS:
#include"unit_tests/containers/lss_matrix_2d_t.hpp"
#include"unit_tests/containers/lss_matrix_3d_t.hpp"
// TESTS FOR SPARSE_SOLVERS:
#include"unit_tests/sparse_solvers/lss_core_cuda_solver_t.hpp"
#include"unit_tests/sparse_solvers/lss_core_sor_solver_cuda_t.hpp"
#include"unit_tests/sparse_solvers/lss_core_sor_solver_t.hpp"
#include"unit_tests/sparse_solvers/lss_cuda_solver_t.hpp"
#include"unit_tests/sparse_solvers/lss_double_sweep_solver_t.hpp"
#include"unit_tests/sparse_solvers/lss_thomas_lu_solver_t.hpp"
#include"unit_tests/sparse_solvers/lss_sor_solver_t.hpp"
#include"unit_tests/sparse_solvers/lss_sor_solver_cuda_t.hpp"
#include"unit_tests/sparse_solvers/lss_karawia_solver_t.hpp"
// TESTS FOR 2ND DEGREE ODE PROBLEMS:
#include"unit_tests/ode_solvers/second_degree/lss_odes_2_degree_t.hpp"
// TESTS FOR 1D PDE PROBLEMS:
#include"unit_tests/pde_solvers/1d/lss_advection_equation_t.hpp"
#include"unit_tests/pde_solvers/1d/lss_black_scholes_equation_t.hpp"
#include"unit_tests/pde_solvers/1d/lss_pure_heat_equation_t.hpp"
#include"unit_tests/pde_solvers/1d/lss_pure_wave_equation_t.hpp"
// TESTS FOR 2D PDE PROBLEMS:
#include"unit_tests/pde_solvers/2d/lss_heston_equation_t.hpp"
#include"unit_tests/pde_solvers/2d/lss_sabr_equation_t.hpp"



int main()
{
  // ======================================================
  // ================== lss_matrix_2d_t ===================
  // ======================================================
  // basic_matrix_2d();
  // slice_matrix_2d();
  // 
  // ======================================================
  // ================== lss_matrix_3d_t ===================
  // ======================================================
  // basic_matrix_3d();
  // slice_matrix_3d();
  // slice_planes_matrix_3d();
  // 
  // ======================================================
  // ============= lss_core_cuda_solver_t =================
  // ======================================================
  // test_device_sparse_qr();
  // test_host_sparse_qr_test();
  // test_dirichlet_bc_bvp_on_host();
  // test_dirichlet_bc_bvp_on_device();
  // test_robin_bc_bvp_on_host();
  // test_robin_bc_bvp_on_device();
  //
  // ======================================================
  // =========== lss_numerical_integration_t ==============
  // ======================================================
  // trapezoid_integral_1_test();
  // 
  // ======================================================
  // =========== lss_numerical_integration_t ==============
  // ======================================================
  // heston_price_test();
  //
  // ======================================================
  // =========== lss_core_sor_solver_cuda_t ===============
  // ======================================================
  // test_device_sor_solver_basic();
  // test_device_sor_solver_ode();
  // 
  // ======================================================
  // ============== lss_core_sor_solver_t =================
  // ======================================================
  // test_host_sor_solver_basic();
  // test_host_sor_solver_ode();
  // 
  // ======================================================
  // ================ lss_cuda_solver_t ===================
  // ======================================================
  // test_device_bvp_dirichlet_bc();
  // test_device_bvp_robin_bc();
  // test_device_bvp_mix_bc();
  // 
  // ======================================================
  // ============ lss_double_sweep_solver_t ===============
  // ======================================================
  // test_dss_bvp_dirichlet_bc();
  // test_dss_bvp_robin_bc();
  // test_dss_bvp_mix_bc();
  // 
  // ======================================================
  // ============== lss_thomas_lu_solver_t ================
  // ======================================================
  // test_thomas_lu_bvp_dirichlet_bc();
  // test_thomas_lu_bvp_robin_bc();
  // test_thomas_lu_bvp_mix_bc();
  // 
  // ======================================================
  // ================= lss_sor_solver_t ===================
  // ======================================================
  // test_host_sor_bvp_dirichlet_bc();
  // test_host_sor_bvp_robin_bc();
  // test_host_sor_bvp_mix_bc();
  // 
  // ======================================================
  // ============== lss_sor_solver_cuda_t =================
  // ======================================================
  // test_device_sor_bvp_dirichlet_bc();
  // test_device_sor_bvp_robin_bc();
  // test_device_sor_bvp_mix_bc();
  // 
  // ======================================================
  // =============== lss_karawia_solver_t =================
  // ======================================================
  // test_karawia_bvp_dirichlet_bc();
  // 
  // ======================================================
  // =============== lss_odes_2_degree_t ==================
  // ======================================================
  // test_impl_simple_ode_dirichlet_bc_cuda_solver_device();
  // test_impl_simple_ode_dirichlet_neumann_bc_cuda_solver_device();
  // test_impl_simple_ode_dirichlet_robin_bc_cuda_solver_device();
  // test_impl_simple_ode_neumann_robin_bc_cuda_solver_device();
  // test_impl_simple_ode1_neumann_robin_bc_cuda_solver_device();
  // 
  // ======================================================
  // =========== lss_advection_equation_t =================
  // ======================================================
  // test_impl_adv_diff_equation_dirichlet_bc_cuda_solver_device_qr();
  // test_impl_adv_diff_equation_dirichlet_bc_sor_solver_device();
  // test_impl_adv_diff_equation_dirichlet_bc_sor_solver_host();
  // test_impl_adv_diff_equation_dirichlet_bc_cuda_solver_host_qr();
  // test_impl_adv_diff_equation_dirichlet_bc_double_sweep_solver();
  // test_impl_adv_diff_equation_dirichlet_bc_thomas_lu_solver();
  // 
  // ======================================================
  // =========== lss_black_scholes_equation_t =============
  // ======================================================
  // test_impl_black_scholes_equation_dirichlet_bc_cuda_solver_device_qr();
  // test_impl_black_scholes_equation_dirichlet_bc_sor_solver_device();
  // test_impl_black_scholes_equation_dirichlet_bc_sor_solver_host();
  // test_impl_black_scholes_equation_dirichlet_bc_double_sweep_solver();
  // test_impl_black_scholes_equation_dirichlet_bc_thomas_lu_solver();
  // test_impl_fwd_black_scholes_equation_dirichlet_bc_cuda_solver_device_qr();
  // test_impl_black_scholes_equation_dirichlet_bc_thomas_lu_solver_stepping();
  // test_expl_black_scholes_equation_dirichlet_bc_ade();
  // 
  // ======================================================
  // ================== lss_pure_heat_equation_t ==========
  // ======================================================
  // test_impl_pure_heat_equation_dirichlet_bc_cuda_solver_device_qr();
  // test_impl_pure_heat_equation_dirichlet_bc_sor_solver_device();
  // test_impl_pure_heat_equation_dirichlet_bc_sor_solver_host();
  // test_impl_pure_heat_equation_dirichlet_bc_cuda_solver_host_qr();
  // test_impl_pure_heat_equation_dirichlet_bc_double_sweep_solver();
  // test_impl_pure_heat_equation_dirichlet_bc_thomas_lu_solver();
  // test_impl_pure_heat_equation_neumann_bc_cuda_solver_device_qr();
  // test_impl_pure_heat_equation_neumann_bc_thomas_lu_solver();
  // test_impl_pure_heat_equation_neumann_bc_double_sweep_solver();
  // test_impl_pure_heat_equation_dirichlet_bc_cuda_solver_device_qr_stepping();
  // test_impl_pure_heat_equation_source_dirichlet_bc_cuda_solver_device();
  // test_impl_pure_heat_equation_source_dirichlet_bc_sor_solver_device();
  // test_expl_pure_heat_equation_dirichlet_bc_ade();
  // test_expl_pure_heat_equation_neumann_bc_euler();
  // test_expl_pure_heat_equation_dirichlet_bc_device();
  // 
  // ======================================================
  // =============== lss_pure_wave_equation_t =============
  // ======================================================
  // test_impl_pure_wave_equation_dirichlet_bc_cuda_solver_device_qr();
  // test_impl_pure_wave_equation_dirichlet_bc_cuda_solver_device_sor();
  // test_impl_pure_wave_equation_dirichlet_bc_cuda_solver_host_sor();
  // test_impl_pure_wave_equation_dirichlet_bc_solver_host_double_sweep();
  // test_impl_pure_wave_equation_dirichlet_bc_solver_host_lu();
  // test_impl_wave_equation_dirichlet_bc_solver_host_lu();
  // test_impl_damped_wave_equation_dirichlet_bc_solver_host_double_sweep();
  // test_impl_pure_wave_equation_neumann_bc_cuda_solver_device_qr();
  // test_expl_pure_wave_equation_dirichlet_bc_cuda_host_solver();
  // test_expl_pure_wave_equation_dirichlet_bc_cuda_device_solver();
  //  
  // ======================================================
  // ================ lss_heston_equation_t ===============
  // ======================================================
  // test_impl_heston_equation_cuda_qr_solver();
  // test_impl_heston_equation_thomas_lu_solver();
  // 
  // ======================================================
  // ================ lss_heston_equation_t ===============
  // ======================================================
  // test_impl_sabr_equation_double_sweep_solver();
  //
  // ======================================================
  // ================ lss_sabr_equation_t =================
  // ======================================================
  // test_impl_sabr_equation_double_sweep_solver();
  // test_impl_sabr_equation_double_sweep_solver_stepping();
  // 
  // ======================================================
  // ===================== lss_xml_t ======================
  // ======================================================
  // test_impl_simple_ode_thomes_lu_qr_xml();
  // test_impl_bs_thomas_lu_crv_xml();
  // test_impl_bs_thomas_lu_srf_xml();
  // test_impl_ph_dirichlet_bvp_device_qr_srf_xml();
  // test_expl_ph_neumann_bvp_euler_srf_xml();
  // test_impl_adv_thomas_lu_srf_xml();      
  // test_impl_pw_dirichlet_bvp_cuda_device_qr_srf_xml();
  // test_impl_w_dirichlet_bvp_host_lu_srf_xml();
  // test_impl_w_bvp_host_dss_srf_xml();
  // test_impl_pw_neumann_bvp_cuda_device_qr_srf_xml();
  // test_impl_heston_cuda_qr_cn_dr_srf_xml();
  // test_impl_heston_thomas_lu_cn_srf_xml();
  // test_impl_sabr_double_sweep_cn_srf_xml();
  // test_impl_heston_thomas_lu_dr_cn_srf_xml();
  // test_impl_heston_dss_cs_cn_srf_xml();
  // test_impl_heston_thomas_lu_mcs_cn_srf_xml();
  // test_impl_heston_thomas_lu_hw_cn_srf_xml();
  // test_expl_heston_host_euler_srf_xml();
  // test_expl_heston_device_euler_srf_xml();
  // test_expl_sabr_host_euler_srf_xml();
  // test_expl_sabr_device_euler_srf_xml();

    std::cin.get();
    std::cin.get();
    return 0;
}


