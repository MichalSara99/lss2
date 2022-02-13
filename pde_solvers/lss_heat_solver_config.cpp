#include "lss_heat_solver_config.hpp"

namespace lss_pde_solvers
{

heat_implicit_solver_config::heat_implicit_solver_config(memory_space_enum const &memory_space,
                                                         traverse_direction_enum const &traverse_direction,
                                                         tridiagonal_method_enum const &tridiagonal_method,
                                                         factorization_enum const &tridiagonal_factorization,
                                                         implicit_pde_scheme_ptr const &pde_scheme_ptr)
    : pde_implicit_solver_config{memory_space, traverse_direction, tridiagonal_method, tridiagonal_factorization}
{
    LSS_ASSERT(pde_scheme_ptr != nullptr, "heat_implicit_solver_config: implicit_pde_scheme must not be empty");
    implicit_pde_scheme_value_ = pde_scheme_ptr->value();
}

heat_implicit_solver_config ::~heat_implicit_solver_config()
{
}

double heat_implicit_solver_config::implicit_pde_scheme_value() const
{
    return implicit_pde_scheme_value_;
}

heat_explicit_solver_config::heat_explicit_solver_config(memory_space_enum const &memory_space,
                                                         traverse_direction_enum const &traverse_direction,
                                                         explicit_pde_schemes_enum const &explicit_pde_scheme)
    : pde_explicit_solver_config{memory_space, traverse_direction}, explicit_pde_scheme_{explicit_pde_scheme}
{
}

heat_explicit_solver_config::~heat_explicit_solver_config()
{
}

explicit_pde_schemes_enum heat_explicit_solver_config::explicit_pde_scheme() const
{
    return explicit_pde_scheme_;
}
} // namespace lss_pde_solvers
