#include "lss_ode_discretization_config.hpp"

namespace lss_ode_solvers
{
ode_discretization_config::ode_discretization_config(range_ptr const &space_range,
                                                     std::size_t const &number_of_space_points)
    : lss_discretization::discretization_config_1d(space_range, number_of_space_points)
{
}

ode_discretization_config ::~ode_discretization_config()
{
}

} // namespace lss_ode_solvers
