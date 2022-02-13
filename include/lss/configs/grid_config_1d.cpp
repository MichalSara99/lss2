#include "grid_config_1d.hpp"

namespace lss
{

grid_config_1d_builder::grid_config_1d_builder()
{
}

grid_config_1d_builder &grid_config_1d_builder::accumulation_point(double accumulation_point)
{
    accumulation_point_ = accumulation_point;
    return *this;
}

grid_config_1d_builder &grid_config_1d_builder::alpha_scale(double alpha_scale)
{
    alpha_scale_ = alpha_scale;
    return *this;
}

grid_config_1d_builder &grid_config_1d_builder::grid_type(grid grid)
{
    grid_ = grid;
    return *this;
}

grid_config_1d_ptr grid_config_1d_builder::build()
{
    return std::make_shared<grid_config_1d>(accumulation_point_, alpha_scale_, grid_);
}

} // namespace lss
