#include "heat_implicit_solver_config.hpp"

namespace lss
{

heat_implicit_solver_config_builder::heat_implicit_solver_config_builder()
{
}

heat_implicit_solver_config_builder &heat_implicit_solver_config_builder::memory(memory_space memory)
{
    memory_space_ = memory;
    return *this;
}

heat_implicit_solver_config_builder &heat_implicit_solver_config_builder::traverse(traverse_direction traverse)
{
    traverse_direction_ = traverse;
    return *this;
}

heat_implicit_solver_config_builder &heat_implicit_solver_config_builder::method(tridiagonal_method method)
{
    tridiagonal_method_ = method;
    return *this;
}

heat_implicit_solver_config_builder &heat_implicit_solver_config_builder::tridiagonal_factorization(
    factorization tridiagonal_factorization)
{
    tridiagonal_factorization_ = tridiagonal_factorization;
    return *this;
}

heat_implicit_solver_config_builder &heat_implicit_solver_config_builder::implicit_pde_scheme(
    implicit_pde_scheme_ptr implicit_pde_scheme)
{
    implicit_pde_scheme_ptr_ = implicit_pde_scheme;
    return *this;
}

heat_implicit_solver_config_ptr heat_implicit_solver_config_builder::build()
{
    return std::make_shared<heat_implicit_solver_config>(memory_space_, traverse_direction_, tridiagonal_method_,
                                                         tridiagonal_factorization_, implicit_pde_scheme_ptr_);
}
} // namespace lss
