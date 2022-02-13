#if !defined(_LSS_HEAT_MODIFIED_CRAIG_SNEYD_METHOD_HPP_)
#define _LSS_HEAT_MODIFIED_CRAIG_SNEYD_METHOD_HPP_

#include <functional>

#include "../../../../boundaries/lss_boundary.hpp"
#include "../../../../common/lss_enumerations.hpp"
#include "../../../../common/lss_macros.hpp"
#include "../../../../containers/lss_container_2d.hpp"
#include "../../../../discretization/lss_discretization.hpp"
#include "../../../../discretization/lss_grid.hpp"
#include "../../../../discretization/lss_grid_config.hpp"
#include "../../../../sparse_solvers/tridiagonal/lss_tridiagonal_solver.hpp"
#include "../../../lss_pde_discretization_config.hpp"
#include "../implicit_coefficients/lss_heston_implicit_coefficients.hpp"
#include "lss_heat_splitting_method.hpp"

namespace lss_pde_solvers
{

namespace two_dimensional
{
using lss_boundary::boundary_2d_pair;
using lss_boundary::boundary_2d_ptr;
using lss_containers::container_2d;
using lss_enumerations::by_enum;
using lss_grids::grid_2d;
using lss_grids::grid_config_2d_ptr;
using lss_utility::container_t;

class implicit_heston_scheme_mcs
{

  public:
    static void rhs_intermed_1(heston_implicit_coefficients_ptr const &cfs, grid_config_2d_ptr const &grid_cfg,
                               std::size_t const &y_index, double const &y, container_2d<by_enum::Row> const &input,
                               double const &time, container_t &solution);

    static void rhs_intermed_1_source(heston_implicit_coefficients_ptr const &cfs, grid_config_2d_ptr const &grid_cfg,
                                      std::size_t const &y_index, double const &y,
                                      container_2d<by_enum::Row> const &input,
                                      container_2d<by_enum::Row> const &inhom_input,
                                      container_2d<by_enum::Row> const &inhom_input_next, double const &time,
                                      container_t &solution);

    static void rhs_intermed_2(heston_implicit_coefficients_ptr const &cfs, grid_config_2d_ptr const &grid_cfg,
                               std::size_t const &x_index, double const &x, container_2d<by_enum::Row> const &input,
                               container_2d<by_enum::Row> const &inhom_input, double const &time,
                               container_t &solution);

    static void rhs_intermed_3(heston_implicit_coefficients_ptr const &cfs, grid_config_2d_ptr const &grid_cfg,
                               std::size_t const &y_index, double const &y, container_2d<by_enum::Row> const &input,
                               container_2d<by_enum::Row> const &inhom_input,
                               container_2d<by_enum::Row> const &inhom_input_next, double const &time,
                               container_t &solution);

    static void rhs(heston_implicit_coefficients_ptr const &cfs, grid_config_2d_ptr const &grid_cfg,
                    std::size_t const &x_index, double const &x, container_2d<by_enum::Row> const &input,
                    container_2d<by_enum::Row> const &inhom_input, double const &time, container_t &solution);
};

/**
    heat_modified_craig_sneyd_method object
 */
class heat_modified_craig_sneyd_method : public heat_splitting_method
{

  private:
    // constant:
    const double cone_ = 1.0;
    // solvers:
    lss_tridiagonal_solver::tridiagonal_solver_ptr solvery_ptr_;
    lss_tridiagonal_solver::tridiagonal_solver_ptr solveru_ptr_;
    // scheme coefficients:
    heston_implicit_coefficients_ptr coefficients_;
    grid_config_2d_ptr grid_cfg_;
    // container:
    container_t low_, diag_, high_, rhs_;

    explicit heat_modified_craig_sneyd_method() = delete;

    void initialize(bool is_heat_source_set);

    void split_0(double const &y, double const &time, container_t &low, container_t &diag, container_t &high);

    void split_1(double const &x, double const &time, container_t &low, container_t &diag, container_t &high);

  public:
    explicit heat_modified_craig_sneyd_method(lss_tridiagonal_solver::tridiagonal_solver_ptr const &solvery_ptr,
                                              lss_tridiagonal_solver::tridiagonal_solver_ptr const &solveru_ptr,
                                              heston_implicit_coefficients_ptr const &coefficients,
                                              grid_config_2d_ptr const &grid_config, bool is_heat_source_set);

    ~heat_modified_craig_sneyd_method();

    heat_modified_craig_sneyd_method(heat_modified_craig_sneyd_method const &) = delete;
    heat_modified_craig_sneyd_method(heat_modified_craig_sneyd_method &&) = delete;
    heat_modified_craig_sneyd_method &operator=(heat_modified_craig_sneyd_method const &) = delete;
    heat_modified_craig_sneyd_method &operator=(heat_modified_craig_sneyd_method &&) = delete;

    void solve(container_2d<by_enum::Row> const &prev_solution, boundary_2d_pair const &horizontal_boundary_pair,
               boundary_2d_pair const &vertical_boundary_pair, double const &time,
               container_2d<by_enum::Row> &solution) override;

    void solve(container_2d<by_enum::Row> const &prev_solution, boundary_2d_pair const &horizontal_boundary_pair,
               boundary_2d_pair const &vertical_boundary_pair, double const &time,
               std::function<double(double, double)> const &heat_source, container_2d<by_enum::Row> &solution) override;
};

} // namespace two_dimensional

} // namespace lss_pde_solvers

#endif ///_LSS_HEAT_MODIFIED_CRAIG_SNEYD_METHOD_HPP_
