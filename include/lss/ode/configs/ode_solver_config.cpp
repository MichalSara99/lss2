#include "ode_solver_config.hpp"

namespace lss
{

ode_implicit_solver_config_builder::ode_implicit_solver_config_builder()
{
}

ode_implicit_solver_config_builder &ode_implicit_solver_config_builder::memory(memory_space memory)
{
    memory_space_ = memory;
    return *this;
}

ode_implicit_solver_config_builder &ode_implicit_solver_config_builder::method(tridiagonal_method method)
{
    tridiagonal_method_ = method;
    return *this;
}

ode_implicit_solver_config_builder &ode_implicit_solver_config_builder::tridiagonal_factorization(
    factorization tridiagonal_factorization)
{
    tridiagonal_factorization_ = tridiagonal_factorization;
    return *this;
}

ode_implicit_solver_config_ptr ode_implicit_solver_config_builder::build()
{
    return std::make_shared<ode_implicit_solver_config>(memory_space_, tridiagonal_method_, tridiagonal_factorization_);
}

} // namespace lss
