#if !defined(_DISCRTEIZATION_CONFIG_2D_HPP_)
#define _DISCRTEIZATION_CONFIG_2D_HPP_

#include "../../../common/lss_enumerations.hpp"
#include "../../../common/lss_macros.hpp"
#include "../../../common/lss_utility.hpp"
#include "../../../pde_solvers/lss_pde_discretization_config.hpp"
#include "../utilities/range.hpp"

namespace lss
{
using discretization_config_2d = lss_pde_solvers::pde_discretization_config_2d;
using discretization_config_2d_ptr = lss_pde_solvers::pde_discretization_config_2d_ptr;

struct discretization_config_2d_builder
{
  private:
    lss_utility::range_ptr space_range_1_;
    lss_utility::range_ptr space_range_2_;
    std::size_t number_of_space_points_1_;
    std::size_t number_of_space_points_2_;
    lss_utility::range_ptr time_range_;
    std::size_t number_of_time_points_;

  public:
    LSS_API explicit discretization_config_2d_builder();

    LSS_API discretization_config_2d_builder &space_range_1(range_ptr const &space_range);

    LSS_API discretization_config_2d_builder &space_range_2(range_ptr const &space_range);

    LSS_API discretization_config_2d_builder &number_of_space_points_1(std::size_t const &number_of_space_points);

    LSS_API discretization_config_2d_builder &number_of_space_points_2(std::size_t const &number_of_space_points);

    LSS_API discretization_config_2d_builder &time_range(range_ptr const &time_range);

    LSS_API discretization_config_2d_builder &number_of_time_points(std::size_t const &number_of_time_points);

    LSS_API discretization_config_2d_ptr build();
};

} // namespace lss
#endif ///_DISCRTEIZATION_CONFIG_2D_HPP_
