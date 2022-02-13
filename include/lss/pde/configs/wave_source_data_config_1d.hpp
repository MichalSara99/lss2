#if !defined(_WAVE_SOURCE_DATA_CONFIG_1D_HPP_)
#define _WAVE_SOURCE_DATA_CONFIG_1D_HPP_

#include <functional>

#include "../../../../common/lss_enumerations.hpp"
#include "../../../../common/lss_macros.hpp"
#include "../../../../common/lss_utility.hpp"
#include "../../../../pde_solvers/lss_wave_data_config.hpp"

namespace lss
{

using wave_source_data_config_1d = lss_pde_solvers::wave_source_data_config_1d;

using wave_source_data_config_1d_ptr = lss_pde_solvers::sptr_t<lss_pde_solvers::wave_source_data_config_1d>;

struct wave_source_data_config_1d_builder
{
  private:
    std::function<double(double, double)> wave_source_;

  public:
    LSS_API explicit wave_source_data_config_1d_builder();

    LSS_API wave_source_data_config_1d_builder &wave_source(std::function<double(double, double)> const &wave_source);

    LSS_API wave_source_data_config_1d_ptr build();
};

} // namespace lss

#endif ///_WAVE_SOURCE_DATA_CONFIG_1D_HPP_
