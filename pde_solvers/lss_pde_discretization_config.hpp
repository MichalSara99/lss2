/**

    @file      lss_pde_discretization_config.hpp
    @brief     PDE discratization configuration
    @details   ~
    @author    Michal Sara
    @date      14.12.2021
    @copyright © Michal Sara, 2021. All right reserved.

**/
#pragma once

#if !defined(_LSS_PDE_DISCRETIZATION_CONFIG_HPP_)
#define _LSS_PDE_DISCRETIZATION_CONFIG_HPP_

#include "../common/lss_enumerations.hpp"
#include "../common/lss_macros.hpp"
#include "../common/lss_range.hpp"
#include "../common/lss_utility.hpp"
#include "../discretization/lss_discretization_config.hpp"

namespace lss_pde_solvers
{

using lss_enumerations::dimension_enum;
using lss_utility::range_ptr;
using lss_utility::sptr_t;

/**
    1D pde_discretization_config structure
 */
struct pde_discretization_config_1d final : public lss_discretization::discretization_config_1d
{
  private:
    range_ptr time_range_;
    std::size_t number_of_time_points_;

    explicit pde_discretization_config_1d() = delete;

  public:
    explicit pde_discretization_config_1d(range_ptr const &space_range, std::size_t const &number_of_space_points,
                                          range_ptr const &time_range, std::size_t const &number_of_time_points);
    ~pde_discretization_config_1d();

    /**
        @brief time range
        @retval
    **/
    LSS_API range_ptr const &time_range() const;

    /**
        @brief number of time points in discretization
        @retval
    **/
    LSS_API std::size_t number_of_time_points() const;

    /**
        @brief value for time step
        @retval
    **/
    LSS_API double time_step() const;
};

/**
    2D pde_discretization_config structure
 */
struct pde_discretization_config_2d
{
  private:
    range_ptr space_range_1_;
    range_ptr space_range_2_;
    range_ptr time_range_;
    std::size_t number_of_space_points_1_;
    std::size_t number_of_space_points_2_;
    std::size_t number_of_time_points_;

    explicit pde_discretization_config_2d() = delete;

  public:
    explicit pde_discretization_config_2d(range_ptr const &space_range_1, range_ptr const &space_range_2,
                                          std::size_t const &number_of_space_points_1,
                                          std::size_t const &number_of_space_points_2, range_ptr const &time_range,
                                          std::size_t const &number_of_time_points);
    ~pde_discretization_config_2d();

    /**
        @brief PDE discretization for first space variable
        @retval
    **/
    LSS_API sptr_t<pde_discretization_config_1d> const pde_discretization_1() const;

    /**
        @brief PDE discretization for second space variable
        @retval
    **/
    LSS_API sptr_t<pde_discretization_config_1d> const pde_discretization_2() const;

    /**
        @brief Space range
        @retval
    **/
    LSS_API std::pair<range_ptr, range_ptr> const space_range() const;

    /**
        @brief time range
        @retval
    **/
    LSS_API range_ptr const &time_range() const;

    /**
        @brief Number of space points in discratization
        @retval
    **/
    LSS_API std::pair<std::size_t, std::size_t> const number_of_space_points() const;

    /**
        @brief Number of time points in discratization
        @retval
    **/
    LSS_API std::size_t number_of_time_points() const;

    /**
        @brief Value for space step
        @retval
    **/
    LSS_API std::pair<double, double> space_step() const;

    /**
        @brief Value for time step
        @retval
    **/
    LSS_API double time_step() const;
};

/**
    3D pde_discretization_config structure
 */
struct pde_discretization_config_3d
{
  private:
    range_ptr space_range_1_;
    range_ptr space_range_2_;
    range_ptr space_range_3_;
    range_ptr time_range_;
    std::size_t number_of_space_points_1_;
    std::size_t number_of_space_points_2_;
    std::size_t number_of_space_points_3_;
    std::size_t number_of_time_points_;

    explicit pde_discretization_config_3d() = delete;

  public:
    explicit pde_discretization_config_3d(range_ptr const &space_range_1, range_ptr const &space_range_2,
                                          range_ptr const &space_range_3, std::size_t const &number_of_space_points_1,
                                          std::size_t const &number_of_space_points_2,
                                          std::size_t const &number_of_space_points_3, range_ptr const &time_range,
                                          std::size_t const &number_of_time_points);
    ~pde_discretization_config_3d();

    /**
        @brief PDE discretization for first space variable
        @retval
    **/
    LSS_API sptr_t<pde_discretization_config_1d> const pde_discretization_1() const;

    /**
        @brief PDE discretization for second space variable
        @retval
    **/
    LSS_API sptr_t<pde_discretization_config_1d> const pde_discretization_2() const;

    /**
        @brief PDE discretization for third space variable
        @retval
    **/
    LSS_API sptr_t<pde_discretization_config_1d> const pde_discretization_3() const;

    /**
        @brief PDE discretization for first and second space variable
        @retval
    **/
    LSS_API sptr_t<pde_discretization_config_2d> const pde_discretization_12() const;

    /**
        @brief PDE discretization for second and first space variable
        @retval
    **/
    LSS_API sptr_t<pde_discretization_config_2d> const pde_discretization_21() const;

    /**
        @brief PDE discretization for first and third space variable
        @retval
    **/
    LSS_API sptr_t<pde_discretization_config_2d> const pde_discretization_13() const;

    /**
        @brief PDE discretization for third and first space variable
        @retval
    **/
    LSS_API sptr_t<pde_discretization_config_2d> const pde_discretization_31() const;

    /**
        @brief PDE discretization for second and third space variable
        @retval
    **/
    LSS_API sptr_t<pde_discretization_config_2d> const pde_discretization_23() const;

    /**
        @brief PDE discretization for third and second space variable
        @retval
    **/
    LSS_API sptr_t<pde_discretization_config_2d> const pde_discretization_32() const;

    /**
        @brief Space ranges
        @retval
    **/
    LSS_API std::tuple<range_ptr, range_ptr, range_ptr> const space_range() const;

    /**
        @brief time range
        @retval
    **/
    LSS_API range_ptr const &time_range() const;

    /**
        @brief  Number of space points in discratiozation
        @retval
    **/
    LSS_API std::tuple<std::size_t, std::size_t, std::size_t> const number_of_space_points() const;

    /**
        @brief Number of time points in discretization
        @retval
    **/
    LSS_API std::size_t number_of_time_points() const;

    /**
        @brief Values for space steps
        @retval
    **/
    LSS_API std::tuple<double, double, double> space_step() const;

    /**
        @brief Value for time step
        @retval
    **/
    LSS_API double time_step() const;
};

using pde_discretization_config_1d_ptr = sptr_t<pde_discretization_config_1d>;
using pde_discretization_config_2d_ptr = sptr_t<pde_discretization_config_2d>;
using pde_discretization_config_3d_ptr = sptr_t<pde_discretization_config_3d>;

} // namespace lss_pde_solvers

#endif ///_LSS_PDE_DISCRETIZATION_CONFIG_HPP_
