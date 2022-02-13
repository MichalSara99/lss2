/**

    @file      lss_robin_boundary.hpp
    @brief     represents general 1D1 2D and 3D Robin boundary
    @details   ~
    @author    Michal Sara
    @date      14.12.2021
    @copyright © Michal Sara, 2021. All right reserved.

**/
#pragma once

#if !defined(_LSS_ROBIN_BOUNDARY_HPP_)
#define _LSS_ROBIN_BOUNDARY_HPP_

#include "../common/lss_macros.hpp"
#include "../common/lss_utility.hpp"
#include "lss_boundary.hpp"

namespace lss_boundary
{

using lss_utility::sptr_t;
/**

    @class   robin_boundary_1d
    @brief   represents general 1D Robin boundary
    @details ~

**/
class robin_boundary_1d final : public boundary_1d
{
  protected:
    robin_boundary_1d() = delete;

  public:
    explicit robin_boundary_1d(double linear_value, double value);

    explicit robin_boundary_1d(const std::function<double(double)> &linear_value,
                               const std::function<double(double)> &value);

    LSS_API double linear_value() const;
    LSS_API double value() const override;

    LSS_API double linear_value(double time) const;
    LSS_API double value(double time) const override;
};
using robin_boundary_1d_ptr = sptr_t<robin_boundary_1d>;

/**

    @class   robin_boundary_2d
    @brief   represents general 2D Robin boundary
    @details ~

**/
class robin_boundary_2d final : public boundary_2d
{
  protected:
    robin_boundary_2d() = delete;

  public:
    explicit robin_boundary_2d(const std::function<double(double, double)> &linear_value,
                               const std::function<double(double, double)> &value);

    LSS_API double linear_value(double time, double space_arg) const;
    LSS_API double value(double time, double space_arg) const override;
};

using robin_boundary_2d_ptr = sptr_t<robin_boundary_2d>;

/**

    @class   robin_boundary_3d
    @brief   represents general 3D Robin boundary
    @details ~

**/
class robin_boundary_3d final : public boundary_3d
{
  protected:
    robin_boundary_3d() = delete;

  public:
    explicit robin_boundary_3d(const std::function<double(double, double, double)> &linear_value,
                               const std::function<double(double, double, double)> &value);

    LSS_API double linear_value(double time, double space_1_arg, double space_2_arg) const;
    LSS_API double value(double time, double space_1_arg, double space_2_arg) const override;
};

using robin_boundary_2d_ptr = sptr_t<robin_boundary_2d>;
} // namespace lss_boundary

#endif ///_LSS_ROBIN_BOUNDARY_1D_HPP_
