#include "lss_heat_explicit_solver_method_2d.hpp"

namespace lss_pde_solvers
{

namespace two_dimensional
{

heat_explicit_solver_method_2d::heat_explicit_solver_method_2d(grid_config_2d_ptr const &grid_config)
    : grid_cfg_{grid_config}
{
}

heat_explicit_solver_method_2d ::~heat_explicit_solver_method_2d()
{
}

} // namespace two_dimensional

} // namespace lss_pde_solvers
