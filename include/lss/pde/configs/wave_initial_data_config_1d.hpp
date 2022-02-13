#if !defined(_WAVE_INITIAL_DATA_CONFIG_1D_HPP_)
#define _WAVE_INITIAL_DATA_CONFIG_1D_HPP_

#include <functional>

#include "../../../../common/lss_enumerations.hpp"
#include "../../../../common/lss_macros.hpp"
#include "../../../../common/lss_utility.hpp"
#include "../../../../pde_solvers/lss_wave_data_config.hpp"

namespace lss
{

using wave_initial_data_config_1d = lss_pde_solvers::wave_initial_data_config_1d;
using wave_initial_data_config_1d_ptr = lss_pde_solvers::sptr_t<lss_pde_solvers::wave_initial_data_config_1d>;

struct wave_initial_data_config_1d_builder
{
  private:
    std::function<double(double)> first_initial_condition_;
    std::function<double(double)> second_initial_condition_;

  public:
    LSS_API explicit wave_initial_data_config_1d_builder();

    LSS_API wave_initial_data_config_1d_builder &first_condition(std::function<double(double)> const &first_condition);

    LSS_API wave_initial_data_config_1d_builder &second_condition(
        std::function<double(double)> const &second_condition);

    LSS_API wave_initial_data_config_1d_ptr build();
};

} // namespace lss

#endif ///_WAVE_INITIAL_DATA_CONFIG_1D_HPP_
