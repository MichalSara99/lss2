#if !defined(_HEAT_DATA_CONFID_2D_HPP_)
#define _HEAT_DATA_CONFID_2D_HPP_

#include <functional>

#include "../../../../common/lss_enumerations.hpp"
#include "../../../../common/lss_macros.hpp"
#include "../../../../common/lss_utility.hpp"
#include "../../../../pde_solvers/lss_heat_data_config.hpp"
#include "heat_coefficient_data_config_2d.hpp"
#include "heat_initial_data_config_2d.hpp"
#include "heat_source_data_config_2d.hpp"

namespace lss
{

using heat_data_config_2d = lss_pde_solvers::heat_data_config_2d;
using heat_data_config_2d_ptr = lss_pde_solvers::sptr_t<lss_pde_solvers::heat_data_config_2d>;

struct heat_data_config_2d_builder
{
  private:
    heat_coefficient_data_config_2d_ptr coefficient_data_config_;
    heat_initial_data_config_2d_ptr initial_data_config_;
    heat_source_data_config_2d_ptr source_data_config_;

  public:
    LSS_API explicit heat_data_config_2d_builder();

    LSS_API heat_data_config_2d_builder &coefficient_data_config(
        heat_coefficient_data_config_2d_ptr const &coefficient_data_config);

    LSS_API heat_data_config_2d_builder &initial_data_config(
        heat_initial_data_config_2d_ptr const &initial_data_config);

    LSS_API heat_data_config_2d_builder &source_data_config(heat_source_data_config_2d_ptr const &source_data_config);

    LSS_API heat_data_config_2d_ptr build();
};
} // namespace lss

#endif //_HEAT_DATA_CONFID_2D_HPP_
