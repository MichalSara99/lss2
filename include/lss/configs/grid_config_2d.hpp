#if !defined(_GRID_CONFIG_2D_HPP_)
#define _GRID_CONFIG_2D_HPP_

#include "../../../common/lss_enumerations.hpp"
#include "../../../common/lss_macros.hpp"
#include "../../../common/lss_utility.hpp"
#include "../../../discretization/lss_grid_config_hints.hpp"

namespace lss
{
using grid_config_2d_ptr = lss_grids::grid_config_hints_2d_ptr;
using grid_config_2d = lss_grids::grid_config_hints_2d;
using grid = lss_enumerations::grid_enum;

struct grid_config_2d_builder
{
  private:
    double accumulation_point_;
    double alpha_scale_;
    double beta_scale_;
    grid grid_;

  public:
    LSS_API explicit grid_config_2d_builder();

    LSS_API grid_config_2d_builder &accumulation_point(double accumulation_point);

    LSS_API grid_config_2d_builder &alpha_scale(double alpha_scale);

    LSS_API grid_config_2d_builder &beta_scale(double beta_scale);

    LSS_API grid_config_2d_builder &grid_type(grid grid);

    LSS_API grid_config_2d_ptr build();
};

} // namespace lss

#endif ///_GRID_CONFIG_2D_HPP_
