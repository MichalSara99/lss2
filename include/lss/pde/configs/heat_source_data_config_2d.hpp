#if !defined(_HEAT_SOURCE_DATA_CONFIG_2D_HPP_)
#define _HEAT_SOURCE_DATA_CONFIG_2D_HPP_

#include <functional>

#include "../../../../common/lss_enumerations.hpp"
#include "../../../../common/lss_macros.hpp"
#include "../../../../common/lss_utility.hpp"
#include "../../../../pde_solvers/lss_heat_data_config.hpp"

namespace lss
{

using heat_source_data_config_2d = lss_pde_solvers::heat_source_data_config_2d;

using heat_source_data_config_2d_ptr = lss_pde_solvers::sptr_t<lss_pde_solvers::heat_source_data_config_2d>;

struct heat_source_data_config_2d_builder
{
  private:
    std::function<double(double, double, double)> heat_source_;

  public:
    LSS_API explicit heat_source_data_config_2d_builder();

    LSS_API heat_source_data_config_2d_builder &heat_source(
        std::function<double(double, double, double)> const &heat_source);

    LSS_API heat_source_data_config_2d_ptr build();
};

} // namespace lss

#endif ///_HEAT_SOURCE_DATA_CONFIG_2D_HPP_
