#include "lss_wave_solver_config.hpp"

namespace lss_pde_solvers
{

wave_implicit_solver_config::wave_implicit_solver_config(memory_space_enum const &memory_space,
                                                         traverse_direction_enum const &traverse_direction,
                                                         tridiagonal_method_enum const &tridiagonal_method,
                                                         factorization_enum const &tridiagonal_factorization)
    : pde_implicit_solver_config{memory_space, traverse_direction, tridiagonal_method, tridiagonal_factorization}
{
}

wave_implicit_solver_config::~wave_implicit_solver_config()
{
}

wave_explicit_solver_config::wave_explicit_solver_config(memory_space_enum const &memory_space,
                                                         traverse_direction_enum const &traverse_direction)
    : pde_explicit_solver_config{memory_space, traverse_direction}
{
}

wave_explicit_solver_config ::~wave_explicit_solver_config()
{
}

} // namespace lss_pde_solvers
