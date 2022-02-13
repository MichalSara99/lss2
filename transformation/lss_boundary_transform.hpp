/**

    @file      lss_boundary_transform.hpp
    @brief     Represents 1D boundary transformations
    @details   ~
    @author    Michal Sara
    @date      14.12.2021
    @copyright © Michal Sara, 2021. All right reserved.

**/
#pragma once
#if !defined(_LSS_BOUNDARY_TRANSFORM_HPP_)
#define _LSS_BOUNDARY_TRANSFORM_HPP_

#include <functional>

#include "../boundaries/lss_boundary.hpp"
#include "../common/lss_utility.hpp"
#include "../discretization/lss_grid_transform_config.hpp"

namespace lss_transformation
{

using lss_boundary::boundary_1d_pair;
using lss_boundary::boundary_1d_ptr;
using lss_grids::grid_transform_config_1d_ptr;
using lss_utility::sptr_t;

/**
   1D boundary_transform structure
 */
struct boundary_transform_1d
{
  private:
    boundary_1d_pair pair_ptr_;

    void initialize(boundary_1d_pair const &boundary_pair, grid_transform_config_1d_ptr const grid_transform_config);

    explicit boundary_transform_1d() = delete;

  public:
    explicit boundary_transform_1d(boundary_1d_pair const &boundary_pair,
                                   grid_transform_config_1d_ptr const grid_transform_config);

    ~boundary_transform_1d();

    boundary_1d_pair const &boundary_pair() const;
};

using boundary_transform_1d_ptr = sptr_t<boundary_transform_1d>;

} // namespace lss_transformation

#endif ///_LSS_BOUNDARY_TRANSFORM_HPP_
