#include "discretization_config_2d.hpp"

namespace lss
{

discretization_config_2d_builder::discretization_config_2d_builder()
{
}

discretization_config_2d_builder &discretization_config_2d_builder::space_range_1(range_ptr const &space_range)
{
    space_range_1_ = space_range;
    return *this;
}

discretization_config_2d_builder &discretization_config_2d_builder::space_range_2(range_ptr const &space_range)
{
    space_range_2_ = space_range;
    return *this;
}

discretization_config_2d_builder &discretization_config_2d_builder::number_of_space_points_1(
    std::size_t const &number_of_space_points)
{
    number_of_space_points_1_ = number_of_space_points;
    return *this;
}

discretization_config_2d_builder &discretization_config_2d_builder::number_of_space_points_2(
    std::size_t const &number_of_space_points)
{
    number_of_space_points_2_ = number_of_space_points;
    return *this;
}

discretization_config_2d_builder &discretization_config_2d_builder::time_range(range_ptr const &time_range)
{
    time_range_ = time_range;
    return *this;
}

discretization_config_2d_builder &discretization_config_2d_builder::number_of_time_points(
    std::size_t const &number_of_time_points)
{
    number_of_time_points_ = number_of_time_points;
    return *this;
}

discretization_config_2d_ptr discretization_config_2d_builder::build()
{
    return std::make_shared<discretization_config_2d>(space_range_1_, space_range_2_, number_of_space_points_1_,
                                                      number_of_space_points_2_, time_range_, number_of_time_points_);
}

} // namespace lss
