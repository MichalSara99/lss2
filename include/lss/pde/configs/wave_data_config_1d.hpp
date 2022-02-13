#if !defined(_WAVE_DATA_CONFID_1D_HPP_)
#define _WAVE_DATA_CONFID_1D_HPP_

#include <functional>

#include "../../../../common/lss_enumerations.hpp"
#include "../../../../common/lss_macros.hpp"
#include "../../../../common/lss_utility.hpp"
#include "../../../../pde_solvers/lss_heat_data_config.hpp"
#include "wave_coefficient_data_config_1d.hpp"
#include "wave_initial_data_config_1d.hpp"
#include "wave_source_data_config_1d.hpp"

namespace lss
{

using wave_data_config_1d = lss_pde_solvers::wave_data_config_1d;
using wave_data_config_1d_ptr = lss_pde_solvers::sptr_t<lss_pde_solvers::wave_data_config_1d>;

struct wave_data_config_1d_builder
{
  private:
    wave_coefficient_data_config_1d_ptr coefficient_data_config_;
    wave_initial_data_config_1d_ptr initial_data_config_;
    wave_source_data_config_1d_ptr source_data_config_;

  public:
    LSS_API explicit wave_data_config_1d_builder();

    LSS_API wave_data_config_1d_builder &coefficient_data_config(
        wave_coefficient_data_config_1d_ptr const &coefficient_data_config);

    LSS_API wave_data_config_1d_builder &initial_data_config(
        wave_initial_data_config_1d_ptr const &initial_data_config);

    LSS_API wave_data_config_1d_builder &source_data_config(wave_source_data_config_1d_ptr const &source_data_config);

    LSS_API wave_data_config_1d_ptr build();
};
} // namespace lss

#endif //_WAVE_DATA_CONFID_1D_HPP_
