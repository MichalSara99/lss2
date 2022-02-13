#if !defined(_LSS_HESTON_BOUNDARY_TRANSFORM_HPP_)
#define _LSS_HESTON_BOUNDARY_TRANSFORM_HPP_

#include <functional>

#include "../../../../boundaries/lss_boundary.hpp"
#include "../../../../common/lss_utility.hpp"
#include "../../../../discretization/lss_grid_transform_config.hpp"

namespace lss_pde_solvers
{

using lss_boundary::boundary_2d_pair;
using lss_boundary::boundary_2d_ptr;
using lss_grids::grid_transform_config_2d_ptr;
using lss_utility::sptr_t;

/**
    heston_boundary_transform structure
 */
struct heston_boundary_transform
{
  private:
    boundary_2d_ptr v_upper_ptr_;
    boundary_2d_pair h_pair_ptr_;

    void initialize(boundary_2d_ptr const &vertical_upper_boundary_ptr,
                    boundary_2d_pair const &horizontal_boundary_pair,
                    grid_transform_config_2d_ptr const grid_transform_config);

    explicit heston_boundary_transform() = delete;

  public:
    explicit heston_boundary_transform(boundary_2d_ptr const &vertical_upper_boundary_ptr,
                                       boundary_2d_pair const &horizontal_boundary_pair,
                                       grid_transform_config_2d_ptr const grid_transform_config);

    ~heston_boundary_transform();

    boundary_2d_ptr const &vertical_upper() const;

    boundary_2d_pair const &horizontal_pair() const;
};

using heston_boundary_transform_ptr = sptr_t<heston_boundary_transform>;

} // namespace lss_pde_solvers

#endif ///_LSS_HESTON_BOUNDARY_TRANSFORM_HPP_
