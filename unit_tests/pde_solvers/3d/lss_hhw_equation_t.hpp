#if !defined(_LSS_HHW_EQUATION_T_HPP_)
#define _LSS_HHW_EQUATION_T_HPP_

// ///////////////////////////////////////////////////////////////////////////
//							HESTON-HULL-WHITE PROBLEMS
// ///////////////////////////////////////////////////////////////////////////

// ===========================================================================
// ========================== IMPLICIT SOLVERS ===============================
// ===========================================================================

extern void impl_hhw_equation_cuda_qr_solver_crank_nicolson();

extern void test_impl_hhw_equation_cuda_qr_solver();

extern void impl_hhw_equation_thomas_lu_solver_crank_nicolson();

extern void impl_hhw_equation_dsssolver_solver_crank_nicolson();

extern void test_impl_hhw_equation_tlu_dss_solver();

#endif //_LSS_HHW_EQUATION_T_HPP_
