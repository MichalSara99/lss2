#if !defined(_LSS_HHW_BOUNDARY_TRANSFORM_HPP_)
#define _LSS_HHW_BOUNDARY_TRANSFORM_HPP_

#include <functional>

#include "../../../../../boundaries/lss_boundary.hpp"
#include "../../../../../common/lss_utility.hpp"
#include "../../../../../discretization/lss_grid_transform_config.hpp"

namespace lss_pde_solvers
{

using lss_boundary::boundary_3d_pair;
using lss_boundary::boundary_3d_ptr;
using lss_grids::grid_transform_config_3d_ptr;
using lss_utility::sptr_t;

/**
    hhw_boundary_transform structure
 */
struct hhw_boundary_transform
{
  private:
    boundary_3d_ptr y_upper_ptr_;
    boundary_3d_pair x_pair_ptr_;
    boundary_3d_pair z_pair_ptr_;

    void initialize(boundary_3d_pair const &x_boundary_pair, boundary_3d_ptr const &y_upper_boundary_ptr,
                    boundary_3d_pair const &z_boundary_pair, grid_transform_config_3d_ptr const grid_transform_config);

    explicit hhw_boundary_transform() = delete;

  public:
    explicit hhw_boundary_transform(boundary_3d_pair const &x_boundary_pair,
                                    boundary_3d_ptr const &y_upper_boundary_ptr,
                                    boundary_3d_pair const &z_boundary_pair,
                                    grid_transform_config_3d_ptr const grid_transform_config);

    ~hhw_boundary_transform();

    boundary_3d_pair const &x_boundary_pair() const;

    boundary_3d_ptr const &y_upper_boundary() const;

    boundary_3d_pair const &z_boundary_pair() const;
};

using hhw_boundary_transform_ptr = sptr_t<hhw_boundary_transform>;

} // namespace lss_pde_solvers

#endif ///_LSS_HHW_BOUNDARY_TRANSFORM_HPP_
