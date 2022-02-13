#if !defined(_DISCRTEIZATION_CONFIG_ODE_HPP_)
#define _DISCRTEIZATION_CONFIG_ODE_HPP_

#include "../../../common/lss_enumerations.hpp"
#include "../../../common/lss_macros.hpp"
#include "../../../common/lss_utility.hpp"
#include "../../../ode_solvers/lss_ode_discretization_config.hpp"
#include "../utilities/range.hpp"

namespace lss
{
using discretization_config_ode = lss_ode_solvers::ode_discretization_config;
using discretization_config_ode_ptr = lss_ode_solvers::ode_discretization_config_ptr;

struct discretization_config_ode_builder
{
  private:
    lss_utility::range_ptr space_range_;
    std::size_t number_of_space_points_;

  public:
    LSS_API explicit discretization_config_ode_builder();

    LSS_API discretization_config_ode_builder &space_range(range_ptr const &space_range);

    LSS_API discretization_config_ode_builder &number_of_space_points(std::size_t const &number_of_space_points);

    LSS_API discretization_config_ode_ptr build();
};

} // namespace lss
#endif ///_DISCRTEIZATION_CONFIG_ODE_HPP_
