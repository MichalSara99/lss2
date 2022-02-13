#if !defined(_DISCRTEIZATION_CONFIG_1D_HPP_)
#define _DISCRTEIZATION_CONFIG_1D_HPP_

#include "../../../common/lss_enumerations.hpp"
#include "../../../common/lss_macros.hpp"
#include "../../../common/lss_utility.hpp"
#include "../../../pde_solvers/lss_pde_discretization_config.hpp"
#include "../utilities/range.hpp"

namespace lss
{
using discretization_config_1d = lss_pde_solvers::pde_discretization_config_1d;
using discretization_config_1d_ptr = lss_pde_solvers::pde_discretization_config_1d_ptr;

struct discretization_config_1d_builder
{
  private:
    lss_utility::range_ptr space_range_;
    std::size_t number_of_space_points_;
    lss_utility::range_ptr time_range_;
    std::size_t number_of_time_points_;

  public:
    LSS_API explicit discretization_config_1d_builder();

    LSS_API discretization_config_1d_builder &space_range(range_ptr const &space_range);

    LSS_API discretization_config_1d_builder &number_of_space_points(std::size_t const &number_of_space_points);

    LSS_API discretization_config_1d_builder &time_range(range_ptr const &time_range);

    LSS_API discretization_config_1d_builder &number_of_time_points(std::size_t const &number_of_time_points);

    LSS_API discretization_config_1d_ptr build();
};

} // namespace lss
#endif ///_DISCRTEIZATION_CONFIG_1D_HPP_
