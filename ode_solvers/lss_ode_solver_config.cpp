#include "lss_ode_solver_config.hpp"

namespace lss_ode_solvers
{

ode_implicit_solver_config::ode_implicit_solver_config(memory_space_enum const &memory_space,
                                                       tridiagonal_method_enum const &tridiagonal_method,
                                                       factorization_enum const &tridiagonal_factorization)
    : memory_space_{memory_space}, tridiagonal_method_{tridiagonal_method}, tridiagonal_factorization_{
                                                                                tridiagonal_factorization}
{
    initialize();
}
ode_implicit_solver_config ::~ode_implicit_solver_config()
{
}

memory_space_enum ode_implicit_solver_config::memory_space() const
{
    return memory_space_;
}

tridiagonal_method_enum ode_implicit_solver_config::tridiagonal_method() const
{
    return tridiagonal_method_;
}

factorization_enum ode_implicit_solver_config::tridiagonal_factorization() const
{
    return tridiagonal_factorization_;
}

void ode_implicit_solver_config::initialize()
{
    if (memory_space_ == memory_space_enum::Device)
    {
        LSS_VERIFY(!(tridiagonal_method_ == tridiagonal_method_enum::DoubleSweepSolver),
                   "No support for Double Sweep Solver on Device");
        LSS_VERIFY(!(tridiagonal_method_ == tridiagonal_method_enum::ThomasLUSolver),
                   "No support for Tomas LU Solver on Device");
    }

    if (tridiagonal_method_ == tridiagonal_method_enum::DoubleSweepSolver)
    {
        LSS_VERIFY(!(tridiagonal_factorization_ == factorization_enum::CholeskyMethod),
                   "No support for Cholesky Method factorization for Double Sweep Solver");
        LSS_VERIFY(!(tridiagonal_factorization_ == factorization_enum::LUMethod),
                   "No support for LU Method factorization for Double Sweep Solver");
        LSS_VERIFY(!(tridiagonal_factorization_ == factorization_enum::QRMethod),
                   "No support for QR Method factorization for Double Sweep "
                   "Solver");
    }

    if (tridiagonal_method_ == tridiagonal_method_enum::ThomasLUSolver)
    {
        LSS_VERIFY(!(tridiagonal_factorization_ == factorization_enum::CholeskyMethod),
                   "No support for Cholesky Method factorization for Thomas "
                   "LU Solver");
        LSS_VERIFY(!(tridiagonal_factorization_ == factorization_enum::LUMethod),
                   "No support for LU Method factorization for Thomas LU Solver");
        LSS_VERIFY(!(tridiagonal_factorization_ == factorization_enum::QRMethod),
                   "No support for QR Method factorization for Thomas LU "
                   "Solver");
    }

    if (tridiagonal_method_ == tridiagonal_method_enum::SORSolver)
    {
        LSS_VERIFY(!(tridiagonal_factorization_ == factorization_enum::CholeskyMethod),
                   "No support for Cholesky Method factorization for SOR"
                   " Solver");
        LSS_VERIFY(!(tridiagonal_factorization_ == factorization_enum::LUMethod),
                   "No support for LU Method factorization for SOR Solver");
        LSS_VERIFY(!(tridiagonal_factorization_ == factorization_enum::QRMethod),
                   "No support for QR Method factorization for SOR "
                   "Solver");
    }
}

} // namespace lss_ode_solvers
