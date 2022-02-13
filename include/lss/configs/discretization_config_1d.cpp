#include "discretization_config_1d.hpp"

namespace lss
{

discretization_config_1d_builder::discretization_config_1d_builder()
{
}

discretization_config_1d_builder &discretization_config_1d_builder::space_range(range_ptr const &space_range)
{
    space_range_ = space_range;
    return *this;
}

discretization_config_1d_builder &discretization_config_1d_builder::number_of_space_points(
    std::size_t const &number_of_space_points)
{
    number_of_space_points_ = number_of_space_points;
    return *this;
}

discretization_config_1d_builder &discretization_config_1d_builder::time_range(range_ptr const &time_range)
{
    time_range_ = time_range;
    return *this;
}

discretization_config_1d_builder &discretization_config_1d_builder::number_of_time_points(
    std::size_t const &number_of_time_points)
{
    number_of_time_points_ = number_of_time_points;
    return *this;
}

discretization_config_1d_ptr discretization_config_1d_builder::build()
{
    return std::make_shared<discretization_config_1d>(space_range_, number_of_space_points_, time_range_,
                                                      number_of_time_points_);
}

} // namespace lss
