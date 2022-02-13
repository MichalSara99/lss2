#include "discretization_config_ode.hpp"

namespace lss
{

discretization_config_ode_builder::discretization_config_ode_builder()
{
}

discretization_config_ode_builder &discretization_config_ode_builder::space_range(range_ptr const &space_range)
{
    space_range_ = space_range;
    return *this;
}

discretization_config_ode_builder &discretization_config_ode_builder::number_of_space_points(
    std::size_t const &number_of_space_points)
{
    number_of_space_points_ = number_of_space_points;
    return *this;
}

discretization_config_ode_ptr discretization_config_ode_builder::build()
{
    return std::make_shared<discretization_config_ode>(space_range_, number_of_space_points_);
}

} // namespace lss
