#if !defined(_GRID_CONFIG_1D_HPP_)
#define _GRID_CONFIG_1D_HPP_

#include "../../../common/lss_enumerations.hpp"
#include "../../../common/lss_macros.hpp"
#include "../../../common/lss_utility.hpp"
#include "../../../discretization/lss_grid_config_hints.hpp"

namespace lss
{
using grid_config_1d_ptr = lss_grids::grid_config_hints_1d_ptr;
using grid_config_1d = lss_grids::grid_config_hints_1d;
using grid = lss_enumerations::grid_enum;

struct grid_config_1d_builder
{
  private:
    double accumulation_point_;
    double alpha_scale_;
    grid grid_;

  public:
    LSS_API explicit grid_config_1d_builder();

    LSS_API grid_config_1d_builder &accumulation_point(double accumulation_point);

    LSS_API grid_config_1d_builder &alpha_scale(double alpha_scale);

    LSS_API grid_config_1d_builder &grid_type(grid grid);

    LSS_API grid_config_1d_ptr build();
};

} // namespace lss

#endif ///_GRID_CONFIG_1D_HPP_
