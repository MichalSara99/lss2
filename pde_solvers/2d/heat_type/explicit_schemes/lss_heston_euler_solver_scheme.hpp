#if !defined(_LSS_HESTON_EULER_SCHEME_HPP_)
#define _LSS_HESTON_EULER_SCHEME_HPP_

#include "../../../../boundaries/lss_boundary.hpp"
#include "../../../../common/lss_enumerations.hpp"
#include "../../../../containers/lss_container_2d.hpp"
#include "../../../../discretization/lss_grid_config.hpp"
#include "../../../lss_pde_discretization_config.hpp"
#include "../boundary_solver/lss_heston_explicit_boundary_solver.hpp"
#include "../explicit_coefficients/lss_heston_euler_coefficients.hpp"
#include "../solver_method/lss_heston_euler_solver_method.hpp"
#include "../time_loop/lss_heston_explicit_time_loop.hpp"

namespace lss_pde_solvers
{

namespace two_dimensional
{

using lss_boundary::boundary_2d_pair;
using lss_boundary::boundary_2d_ptr;
using lss_containers::container_2d;
using lss_enumerations::by_enum;
using lss_enumerations::traverse_direction_enum;
using lss_grids::grid_config_2d_ptr;

class heston_euler_scheme
{
  private:
    heston_euler_coefficients_ptr euler_coeffs_;
    heston_explicit_boundary_solver_ptr heston_boundary_;
    boundary_2d_ptr boundary_ver_;
    boundary_2d_pair boundary_pair_hor_;
    pde_discretization_config_2d_ptr discretization_cfg_;
    grid_config_2d_ptr grid_cfg_;

    bool is_stable(heston_implicit_coefficients_ptr const &coefficients);

    void initialize(heston_implicit_coefficients_ptr const &coefficients);

    explicit heston_euler_scheme() = delete;

  public:
    heston_euler_scheme(heston_implicit_coefficients_ptr const &coefficients,
                        boundary_2d_ptr const &vertical_upper_boundary_ptr,
                        boundary_2d_pair const &horizontal_boundary_pair,
                        pde_discretization_config_2d_ptr const &discretization_config,
                        grid_config_2d_ptr const &grid_config);

    ~heston_euler_scheme();

    void operator()(container_2d<by_enum::Row> &prev_solution, container_2d<by_enum::Row> &next_solution,
                    bool is_heat_sourse_set, std::function<double(double, double, double)> const &heat_source,
                    traverse_direction_enum traverse_dir);
};
} // namespace two_dimensional
} // namespace lss_pde_solvers
#endif ///_LSS_HESTON_EULER_SCHEME_HPP_
