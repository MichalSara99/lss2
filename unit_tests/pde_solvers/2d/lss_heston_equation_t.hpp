#if !defined(_LSS_HESTON_EQUATION_T_HPP_)
#define _LSS_HESTON_EQUATION_T_HPP_

// ///////////////////////////////////////////////////////////////////////////
//							HESTON PROBLEMS
// ///////////////////////////////////////////////////////////////////////////

// ===========================================================================
// ========================== IMPLICIT SOLVERS ===============================
// ===========================================================================

extern void impl_heston_equation_cuda_qr_solver_crank_nicolson();

extern void test_impl_heston_equation_cuda_qr_solver();

extern void impl_heston_equation_thomas_lu_solver_crank_nicolson();

extern void test_impl_heston_equation_thomas_lu_solver();

#endif //_LSS_HESTON_EQUATION_T_HPP_
