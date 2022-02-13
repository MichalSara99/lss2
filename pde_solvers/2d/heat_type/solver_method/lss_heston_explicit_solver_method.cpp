#include "lss_heston_explicit_solver_method.hpp"

namespace lss_pde_solvers
{

namespace two_dimensional
{

heston_explicit_solver_method::heston_explicit_solver_method(grid_config_2d_ptr const &grid_config)
    : grid_cfg_{grid_config}
{
}

heston_explicit_solver_method ::~heston_explicit_solver_method()
{
}

} // namespace two_dimensional

} // namespace lss_pde_solvers
