#include "lss_discretization_config.hpp"

namespace lss_discretization
{

discretization_config_1d::discretization_config_1d(range_ptr const &space_range,
                                                   std::size_t const &number_of_space_points)
    : space_range_{space_range}, number_of_space_points_{number_of_space_points}
{
}

discretization_config_1d ::~discretization_config_1d()
{
}

range_ptr const &discretization_config_1d::space_range() const
{
    return space_range_;
}

std::size_t discretization_config_1d::number_of_space_points() const
{
    return number_of_space_points_;
}

double discretization_config_1d::space_step() const
{
    return ((space_range_->spread()) / static_cast<double>(number_of_space_points_ - 1));
}

} // namespace lss_discretization
