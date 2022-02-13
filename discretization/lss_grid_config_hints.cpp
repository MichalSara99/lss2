#include "lss_grid_config_hints.hpp"

namespace lss_grids
{

grid_config_hints_1d::grid_config_hints_1d(double accumulation_point, double alpha_scale, grid_enum grid_type)
    : accumulation_point_{accumulation_point}, alpha_scale_{alpha_scale}, grid_{grid_type}
{
}

double grid_config_hints_1d::accumulation_point() const
{
    return accumulation_point_;
}

double grid_config_hints_1d::alpha_scale() const
{
    return alpha_scale_;
}

grid_enum grid_config_hints_1d::grid() const
{
    return grid_;
}

grid_config_hints_2d::grid_config_hints_2d(double accumulation_point, double alpha_scale, double beta_scale,
                                           grid_enum grid_type)
    : accumulation_point_{accumulation_point}, alpha_scale_{alpha_scale}, beta_scale_{beta_scale}, grid_{grid_type}
{
}

double grid_config_hints_2d::accumulation_point() const
{
    return accumulation_point_;
}

double grid_config_hints_2d::alpha_scale() const
{
    return alpha_scale_;
}

double grid_config_hints_2d::beta_scale() const
{
    return beta_scale_;
}

grid_enum grid_config_hints_2d::grid() const
{
    return grid_;
}

grid_config_hints_3d::grid_config_hints_3d(double accumulation_point, double alpha_scale, double beta_scale_1,
                                           double beta_scale_2, grid_enum grid_type)
    : accumulation_point_{accumulation_point}, alpha_scale_{alpha_scale}, beta_scale_1_{beta_scale_1},
      beta_scale_2_{beta_scale_2}, grid_{grid_type}
{
}

double grid_config_hints_3d::accumulation_point() const
{
    return accumulation_point_;
}

double grid_config_hints_3d::alpha_scale() const
{
    return alpha_scale_;
}

std::pair<double, double> grid_config_hints_3d::beta_scales() const
{
    return std::make_pair(beta_scale_1_, beta_scale_2_);
}

grid_enum grid_config_hints_3d::grid() const
{
    return grid_;
}

} // namespace lss_grids
