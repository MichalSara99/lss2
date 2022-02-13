#if !defined(_LSS_HESTON_EXPLICIT_BOUNDARY_SOLVER_HPP_)
#define _LSS_HESTON_EXPLICIT_BOUNDARY_SOLVER_HPP_

#include <vector>

#include "../../../../boundaries/lss_boundary.hpp"
#include "../../../../common/lss_enumerations.hpp"
#include "../../../../common/lss_utility.hpp"
#include "../../../../containers/lss_container_2d.hpp"
#include "../../../../discretization/lss_grid_config.hpp"
#include "../implicit_coefficients/lss_heston_implicit_coefficients.hpp"

namespace lss_pde_solvers
{
namespace two_dimensional
{

using lss_boundary::boundary_2d_pair;
using lss_boundary::boundary_2d_ptr;
using lss_containers::container_2d;
using lss_enumerations::by_enum;
using lss_grids::grid_config_2d_ptr;
using lss_utility::container_t;

/**
    explicit_heston_boundary_scheme object
 */
class explicit_heston_boundary_scheme
{

  public:
    static void rhs(heston_implicit_coefficients_ptr const &cfg, grid_config_2d_ptr const &grid_cfg,
                    std::size_t const &y_index, double const &y, boundary_2d_pair const &horizontal_boundary_pair,
                    container_2d<by_enum::Row> const &input, double const &time, container_t &solution);

    static void rhs_source(heston_implicit_coefficients_ptr const &cfg, grid_config_2d_ptr const &grid_cfg,
                           std::size_t const &y_index, double const &y,
                           boundary_2d_pair const &horizontal_boundary_pair, container_2d<by_enum::Row> const &input,
                           container_2d<by_enum::Row> const &inhom_input, double const &time, container_t &solution);
};

/**
    heston_explicit_boundary_solver object
 */
class heston_explicit_boundary_solver
{

  private:
    heston_implicit_coefficients_ptr coefficients_;
    grid_config_2d_ptr grid_cfg_;

    explicit heston_explicit_boundary_solver() = delete;

  public:
    explicit heston_explicit_boundary_solver(heston_implicit_coefficients_ptr const &coefficients,
                                             grid_config_2d_ptr const &grid_config);

    ~heston_explicit_boundary_solver();

    heston_explicit_boundary_solver(heston_explicit_boundary_solver const &) = delete;
    heston_explicit_boundary_solver(heston_explicit_boundary_solver &&) = delete;
    heston_explicit_boundary_solver &operator=(heston_explicit_boundary_solver const &) = delete;
    heston_explicit_boundary_solver &operator=(heston_explicit_boundary_solver &&) = delete;

    void solve(container_2d<by_enum::Row> const &prev_solution, boundary_2d_pair const &horizonatal_boundary_pair,
               boundary_2d_ptr const &vertical_upper_boundary_ptr, double const &time,
               container_2d<by_enum::Row> &solution);

    void solve(container_2d<by_enum::Row> const &prev_solution, boundary_2d_pair const &horizonatal_boundary_pair,
               double const &time, container_2d<by_enum::Row> &solution);
};

using heston_explicit_boundary_solver_ptr = sptr_t<heston_explicit_boundary_solver>;

} // namespace two_dimensional

} // namespace lss_pde_solvers

#endif ///_LSS_HESTON_EXPLICIT_BOUNDARY_SOLVER_HPP_
