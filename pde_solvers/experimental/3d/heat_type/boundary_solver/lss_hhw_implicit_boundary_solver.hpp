#if !defined(_LSS_HHW_IMPLICIT_BOUNDARY_SOLVER_HPP_)
#define _LSS_HHW_IMPLICIT_BOUNDARY_SOLVER_HPP_

#include <vector>

#include "../../../../boundaries/lss_boundary.hpp"
#include "../../../../common/lss_enumerations.hpp"
#include "../../../../common/lss_utility.hpp"
#include "../../../../containers/lss_container_3d.hpp"
#include "../../../../discretization/lss_grid_config.hpp"
#include "../../../../sparse_solvers/tridiagonal/lss_tridiagonal_solver.hpp"
#include "../implicit_coefficients/lss_hhw_implicit_coefficients.hpp"

namespace lss_pde_solvers
{
namespace three_dimensional
{

using lss_boundary::boundary_3d_pair;
using lss_boundary::boundary_3d_ptr;
using lss_containers::container_3d;
using lss_enumerations::by_enum;
using lss_grids::grid_config_3d_ptr;
using lss_utility::container_t;

/**
    implicit_hhw_boundary_scheme object
 */
class implicit_hhw_boundary_scheme
{

  public:
    static void rhs(hhw_implicit_coefficients_ptr const &cfg, grid_config_3d_ptr const &grid_cfg,
                    std::size_t const &y_index, double const &y, boundary_3d_pair const &x_boundary_pair,
                    boundary_3d_pair const &z_boundary_pair, container_3d<by_enum::RowPlane> const &input,
                    double const &time, container_2d<by_enum::Row> &solution);

    static void rhs_source(hhw_implicit_coefficients_ptr const &cfg, grid_config_3d_ptr const &grid_cfg,
                           std::size_t const &y_index, double const &y, boundary_3d_pair const &x_boundary_pair,
                           boundary_3d_pair const &z_boundary_pair, container_3d<by_enum::RowPlane> const &input,
                           container_3d<by_enum::RowPlane> const &inhom_input, double const &time,
                           container_2d<by_enum::Row> &solution);
};

/**
    hhw_implicit_boundary_solver object
 */
class hhw_implicit_boundary_solver
{

  private:
    // constants:
    const double cone_ = 1.0;
    const double ctwo_ = 2.0;
    // solvers:
    lss_tridiagonal_solver::tridiagonal_solver_ptr solveru_ptr_;
    // scheme coefficients:
    hhw_implicit_coefficients_ptr coefficients_;
    grid_config_3d_ptr grid_cfg_;
    // containers:
    container_t low_, diag_, high_, rhs_;

    explicit hhw_implicit_boundary_solver() = delete;

  public:
    explicit hhw_implicit_boundary_solver(lss_tridiagonal_solver::tridiagonal_solver_ptr const &solveru,
                                          hhw_implicit_coefficients_ptr const &coefficients,
                                          grid_config_3d_ptr const &grid_config);

    ~hhw_implicit_boundary_solver();

    hhw_implicit_boundary_solver(hhw_implicit_boundary_solver const &) = delete;
    hhw_implicit_boundary_solver(hhw_implicit_boundary_solver &&) = delete;
    hhw_implicit_boundary_solver &operator=(hhw_implicit_boundary_solver const &) = delete;
    hhw_implicit_boundary_solver &operator=(hhw_implicit_boundary_solver &&) = delete;

    void split(double const &x, double const &y, double const &time, container_t &low, container_t &diag,
               container_t &high);

    void solve(container_3d<by_enum::RowPlane> const &prev_solution, boundary_3d_pair const &x_boundary_pair,
               boundary_3d_ptr const &y_upper_boundary_ptr, boundary_3d_pair const &z_boundary_pair, double const &time,
               container_3d<by_enum::RowPlane> &solution);

    void solve(container_3d<by_enum::RowPlane> const &prev_solution, boundary_3d_pair const &x_boundary_pair,
               boundary_3d_pair const &z_boundary_pair, double const &time, container_3d<by_enum::RowPlane> &solution);
};

using hhw_implicit_boundary_solver_ptr = sptr_t<hhw_implicit_boundary_solver>;

} // namespace three_dimensional

} // namespace lss_pde_solvers

#endif /// _LSS_HHW_IMPLICIT_BOUNDARY_SOLVER_HPP_
