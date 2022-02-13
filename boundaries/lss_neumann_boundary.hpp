/**

    @file      lss_neumann_boundary.hpp
    @brief     represents general 1D, 2D and 3D Neumann boundary
    @details   ~
    @author    Michal Sara
    @date      14.12.2021
    @copyright © Michal Sara, 2021. All right reserved.

**/
#pragma once
#if !defined(_LSS_NEUMANN_BOUNDARY_HPP_)
#define _LSS_NEUMANN_BOUNDARY_HPP_

#include "../common/lss_macros.hpp"
#include "../common/lss_utility.hpp"
#include "lss_boundary.hpp"

namespace lss_boundary
{

using lss_utility::sptr_t;

/**

    @class   neumann_boundary_1d
    @brief   represents general 1D Neumann boundary
    @details ~

**/
class neumann_boundary_1d final : public boundary_1d
{
  protected:
    neumann_boundary_1d() = delete;

  public:
    explicit neumann_boundary_1d(double value);

    explicit neumann_boundary_1d(const std::function<double(double)> &value);

    LSS_API double value(double time) const override;

    LSS_API double value() const override;
};
using neumann_boundary_1d_ptr = sptr_t<neumann_boundary_1d>;

/**

    @class   neumann_boundary_2d
    @brief   represents general 2D Neumann boundary
    @details ~

**/
class neumann_boundary_2d final : public boundary_2d
{
  protected:
    neumann_boundary_2d() = delete;

  public:
    explicit neumann_boundary_2d(const std::function<double(double, double)> &value);

    LSS_API double value(double time, double space_arg) const override;
};

using neumann_boundary_2d_ptr = sptr_t<neumann_boundary_2d>;

/**

    @class   neumann_boundary_3d
    @brief   represents general 3D Neumann boundary
    @details ~

**/
class neumann_boundary_3d final : public boundary_3d
{
  protected:
    neumann_boundary_3d() = delete;

  public:
    explicit neumann_boundary_3d(const std::function<double(double, double, double)> &value);

    LSS_API double value(double time, double space_1_arg, double space_2_arg) const override;
};

using neumann_boundary_3d_ptr = sptr_t<neumann_boundary_3d>;

} // namespace lss_boundary

#endif ///_LSS_NEUMANN_BOUNDARY_HPP_
