#if !defined(_HEAT_DATA_CONFID_1D_HPP_)
#define _HEAT_DATA_CONFID_1D_HPP_

#include <functional>

#include "../../../../common/lss_enumerations.hpp"
#include "../../../../common/lss_macros.hpp"
#include "../../../../common/lss_utility.hpp"
#include "../../../../pde_solvers/lss_heat_data_config.hpp"
#include "heat_coefficient_data_config_1d.hpp"
#include "heat_initial_data_config_1d.hpp"
#include "heat_source_data_config_1d.hpp"

namespace lss
{

using heat_data_config_1d = lss_pde_solvers::heat_data_config_1d;
using heat_data_config_1d_ptr = lss_pde_solvers::sptr_t<lss_pde_solvers::heat_data_config_1d>;

struct heat_data_config_1d_builder
{
  private:
    heat_coefficient_data_config_1d_ptr coefficient_data_config_;
    heat_initial_data_config_1d_ptr initial_data_config_;
    heat_source_data_config_1d_ptr source_data_config_;

  public:
    LSS_API explicit heat_data_config_1d_builder();

    LSS_API heat_data_config_1d_builder &coefficient_data_config(
        heat_coefficient_data_config_1d_ptr const &coefficient_data_config);

    LSS_API heat_data_config_1d_builder &initial_data_config(
        heat_initial_data_config_1d_ptr const &initial_data_config);

    LSS_API heat_data_config_1d_builder &source_data_config(heat_source_data_config_1d_ptr const &source_data_config);

    LSS_API heat_data_config_1d_ptr build();
};
} // namespace lss

#endif //_HEAT_DATA_CONFID_1D_HPP_
