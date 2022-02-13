#if !defined(_LSS_HEAT_DOUGLAS_RACHFORD_METHOD_3D_HPP_)
#define _LSS_HEAT_DOUGLAS_RACHFORD_METHOD_3D_HPP_

#include <functional>
#include <map>

#include "../../../../../boundaries/lss_boundary.hpp"
#include "../../../../../common/lss_enumerations.hpp"
#include "../../../../../containers/lss_container_2d.hpp"
#include "../../../../../discretization/lss_grid_config.hpp"
#include "../../../../../sparse_solvers/tridiagonal/lss_tridiagonal_solver.hpp"
#include "../implicit_coefficients/lss_hhw_implicit_coefficients.hpp"
#include "lss_heat_splitting_method_3d.hpp"

namespace lss_pde_solvers
{

namespace three_dimensional
{
using lss_boundary::boundary_2d_pair;
using lss_boundary::boundary_2d_ptr;
using lss_containers::container_2d;
using lss_enumerations::by_enum;
using lss_grids::grid_config_2d_ptr;
using lss_utility::container_t;

class implicit_hhw_scheme
{

  public:
    static void rhs_intermed_1(hhw_implicit_coefficients_ptr const &cfs, grid_config_3d_ptr const &grid_cfg,
                               std::size_t const &y_index, double const &y, std::size_t const &z_index, double const &z,
                               container_3d<by_enum::RowPlane> const &input, double const &time, container_t &solution);

    static void rhs_intermed_1_source(hhw_implicit_coefficients_ptr const &cfs, grid_config_3d_ptr const &grid_cfg,
                                      std::size_t const &y_index, double const &y, std::size_t const &z_index,
                                      double const &z, container_3d<by_enum::RowPlane> const &input,
                                      container_3d<by_enum::RowPlane> const &inhom_input,
                                      container_3d<by_enum::RowPlane> const &inhom_input_next, double const &time,
                                      container_t &solution);

    static void rhs_intermed_2(hhw_implicit_coefficients_ptr const &cfs, grid_config_3d_ptr const &grid_cfg,
                               std::size_t const &x_index, double const &x, std::size_t const &z_index, double const &z,
                               container_3d<by_enum::RowPlane> const &input,
                               container_3d<by_enum::RowPlane> const &inhom_input, double const &time,
                               container_t &solution);

    static void rhs(hhw_implicit_coefficients_ptr const &cfs, grid_config_3d_ptr const &grid_cfg,
                    std::size_t const &x_index, double const &x, std::size_t const &y_index, double const &y,
                    container_3d<by_enum::RowPlane> const &input, container_3d<by_enum::RowPlane> const &inhom_input,
                    double const &time, container_t &solution);
};

/**
    heat_douglas_rachford_method_3d object
 */
class heat_douglas_rachford_method_3d : public heat_splitting_method_3d
{

  private:
    // constants:
    const double cone_ = 1.0;
    const double ctwo_ = 2.0;
    // solvers:
    lss_tridiagonal_solver::tridiagonal_solver_ptr solvery1_ptr_;
    lss_tridiagonal_solver::tridiagonal_solver_ptr solvery2_ptr_;
    lss_tridiagonal_solver::tridiagonal_solver_ptr solveru_ptr_;
    // scheme coefficients:
    hhw_implicit_coefficients_ptr coefficients_;
    grid_config_3d_ptr grid_cfg_;
    // containers:
    container_t low_, diag_, high_, rhs_;

    explicit heat_douglas_rachford_method_3d() = delete;

    void initialize(bool is_heat_source_set);

    void split_0(double const &y, double const &z, double const &time, container_t &low, container_t &diag,
                 container_t &high);

    void split_1(double const &x, double const &z, double const &time, container_t &low, container_t &diag,
                 container_t &high);

    void split_2(double const &x, double const &y, double const &time, container_t &low, container_t &diag,
                 container_t &high);

  public:
    explicit heat_douglas_rachford_method_3d(lss_tridiagonal_solver::tridiagonal_solver_ptr const &solvery1_ptr,
                                             lss_tridiagonal_solver::tridiagonal_solver_ptr const &solvery2_ptr,
                                             lss_tridiagonal_solver::tridiagonal_solver_ptr const &solveru_ptr,
                                             hhw_implicit_coefficients_ptr const &coefficients,
                                             grid_config_3d_ptr const &grid_config, bool is_heat_source_set);

    ~heat_douglas_rachford_method_3d();

    heat_douglas_rachford_method_3d(heat_douglas_rachford_method_3d const &) = delete;
    heat_douglas_rachford_method_3d(heat_douglas_rachford_method_3d &&) = delete;
    heat_douglas_rachford_method_3d &operator=(heat_douglas_rachford_method_3d const &) = delete;
    heat_douglas_rachford_method_3d &operator=(heat_douglas_rachford_method_3d &&) = delete;

    void solve(container_3d<by_enum::RowPlane> const &prev_solution, boundary_3d_pair const &x_boundary_pair,
               boundary_3d_pair const &y_boundary_pair, boundary_3d_pair const &z_boundary_pair, double const &time,
               container_3d<by_enum::RowPlane> &solution) override;

    void solve(container_3d<by_enum::RowPlane> const &prev_solution, boundary_3d_pair const &x_boundary_pair,
               boundary_3d_pair const &y_boundary_pair, boundary_3d_pair const &z_boundary_pair, double const &time,
               std::function<double(double, double, double)> const &heat_source,
               container_3d<by_enum::RowPlane> &solution) override;
};

} // namespace three_dimensional

} // namespace lss_pde_solvers

#endif ///_LSS_HEAT_DOUGLAS_RACHFORD_METHOD_3D_HPP_
