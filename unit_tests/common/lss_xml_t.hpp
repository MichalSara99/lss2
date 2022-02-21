#if !defined(_LSS_XML_T_HPP_)
#define _LSS_XML_T_HPP_

// ODEs

extern void impl_simple_ode_thomes_lu_qr_xml();

extern void test_impl_simple_ode_thomes_lu_qr_xml();

// PDEs
extern void impl_bs_thomas_lu_solver_euler_crv_xml();

extern void impl_bs_thomas_lu_solver_cn_crv_xml();

extern void test_impl_bs_thomas_lu_crv_xml();

extern void impl_bs_thomas_lu_euler_srf_xml();

extern void impl_bs_thomas_lu_cn_srf_xml();

extern void test_impl_bs_thomas_lu_srf_xml();

extern void impl_ph_dirichlet_bvp_device_qr_euler_srf_xml();

extern void impl_ph_dirichlet_bvp_device_qr_cn_srf_xml();

extern void test_impl_ph_dirichlet_bvp_device_qr_srf_xml();

extern void expl_ph_neumann_neumann_bvp_euler_srf_xml();

extern void test_expl_ph_neumann_bvp_euler_srf_xml();

extern void impl_adv_thomas_lu_solver_cn_srf_xml();

extern void test_impl_adv_thomas_lu_srf_xml();

extern void impl_pw_dirichlet_bvp_cuda_device_qr_srf_xml();

extern void test_impl_pw_dirichlet_bvp_cuda_device_qr_srf_xml();

extern void impl_w_dirichlet_bvp_host_lu_srf_xml();

extern void test_impl_w_dirichlet_bvp_host_lu_srf_xml();

extern void impl_w_bvp_host_double_sweep_srf_xml();

extern void test_impl_w_bvp_host_dss_srf_xml();

extern void impl_pw_neumann_bvp_cuda_device_qr_srf_xml();

extern void test_impl_pw_neumann_bvp_cuda_device_qr_srf_xml();

extern void expl_pw_dirichlet_bvp_cuda_host_srf_xml();

extern void test_expl_pw_dirichlet_bvp_cuda_host_host_srf_xml();

extern void impl_heston_cuda_qr_crank_nicolson_dr_srf_xml();

extern void test_impl_heston_cuda_qr_cn_dr_srf_xml();

extern void impl_heston_upoutcall_barrier_thomas_lu_cn_srf_xml_stepping();

extern void test_impl_heston_thomas_lu_cn_srf_xml_stepping();

extern void impl_heston_thomas_lu_cn_srf_xml();

extern void impl_heston_upoutcall_barrier_thomas_lu_cn_srf_xml();

extern void impl_heston_downoutcall_barrier_put_thomas_lu_cn_srf_xml();

extern void impl_heston_upoutput_barrier_put_thomas_lu_cn_srf_xml();

extern void impl_heston_downoutput_barrier_put_thomas_lu_cn_srf_xml();

extern void test_impl_heston_thomas_lu_cn_srf_xml();

extern void impl_sabr_double_sweep_cn_srf_xml();

extern void impl_sabr_double_sweep_cn_srf_xml_stepping();

extern void impl_sabr_barrier_double_sweep_cn_srf_xml();

extern void test_impl_sabr_double_sweep_cn_srf_xml();

extern void impl_heston_thomas_lu_dr_cn_srf_xml();

extern void test_impl_heston_thomas_lu_dr_cn_srf_xml();

extern void impl_heston_dss_solver_cs_cn_srf_xml();

extern void impl_heston_put_dss_solver_cs_cn_srf_xml();

extern void test_impl_heston_dss_cs_cn_srf_xml();

extern void impl_heston_thomas_lu_mcs_cn_srf_xml();

extern void test_impl_heston_thomas_lu_mcs_cn_srf_xml();

extern void impl_heston_thomas_lu_hw_cn_srf_xml();

extern void test_impl_heston_thomas_lu_hw_cn_srf_xml();

// extern void impl_hhw_dsssolver_dr_cn_srf_xml();
//
// extern void impl_hhw_dsssolver_dr_0_66_cn_srf_xml();
//
// extern void test_impl_hhw_dsssolver_dr_cn_srf_xml();

////////////////////////////////////////////////////////////////////////////////////
///
///                             EXPLICIT SOLVERS
///
////////////////////////////////////////////////////////////////////////////////////

extern void expl_heston_host_euler_srf_xml();

extern void test_expl_heston_host_euler_srf_xml();

extern void expl_heston_device_euler_srf_xml();

extern void test_expl_heston_device_euler_srf_xml();

extern void expl_sabr_host_euler_srf_xml();

extern void test_expl_sabr_host_euler_srf_xml();

extern void expl_sabr_device_euler_srf_xml();

extern void test_expl_sabr_device_euler_srf_xml();

#endif ///_LSS_XML_T_HPP_
