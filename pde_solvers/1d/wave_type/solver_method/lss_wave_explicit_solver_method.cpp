#include "lss_wave_explicit_solver_method.hpp"

namespace lss_pde_solvers
{

namespace one_dimensional
{

wave_explicit_solver_method::wave_explicit_solver_method(grid_config_1d_ptr const &grid_config) : grid_cfg_{grid_config}
{
}

wave_explicit_solver_method ::~wave_explicit_solver_method()
{
}

} // namespace one_dimensional

} // namespace lss_pde_solvers
