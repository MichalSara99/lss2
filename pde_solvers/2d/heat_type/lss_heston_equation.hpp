#if !defined(_LSS_HESTON_EQUATION_HPP_)
#define _LSS_HESTON_EQUATION_HPP_

#include <functional>
#include <map>

#include "../../../boundaries/lss_boundary.hpp"
#include "../../../common/lss_enumerations.hpp"
#include "../../../common/lss_macros.hpp"
#include "../../../containers/lss_container_2d.hpp"
#include "../../../containers/lss_container_3d.hpp"
#include "../../../discretization/lss_discretization.hpp"
#include "../../../discretization/lss_grid_config_hints.hpp"
#include "../../../discretization/lss_grid_transform_config.hpp"
#include "../../lss_heat_data_config.hpp"
#include "../../lss_heat_solver_config.hpp"
#include "../../lss_pde_discretization_config.hpp"
#include "../../lss_splitting_method_config.hpp"
#include "../../transformation/lss_heat_data_transform.hpp"
#include "transformation/lss_heston_boundary_transform.hpp"

namespace lss_pde_solvers
{

namespace two_dimensional
{
using lss_boundary::boundary_2d_pair;
using lss_boundary::boundary_2d_ptr;
using lss_containers::container_2d;
using lss_containers::container_3d;
using lss_enumerations::by_enum;
using lss_enumerations::grid_enum;
using lss_grids::grid_config_hints_2d_ptr;
using lss_grids::grid_transform_config_2d_ptr;

using d_2d = discretization_2d<std::vector, std::allocator<double>>;

namespace implicit_solvers
{

/*!
============================================================================
Represents general variable coefficient Heston type equation

u_t = a(t,x,y)*u_xx + b(t,x,y)*u_yy + c(t,x,y)*u_xy + d(t,x,y)*u_x + e(t,x,y)*u_y +
        f(t,x,y)*u + F(t,x,y)

t > 0, x_1 < x < x_2, y_1 < y < y_2

with initial condition:

u(0,x,y) = G(x,y)

or terminal condition:

u(T,x,y) = G(x,y)

horizontal_boundary_pair = S = (S_1,S_2) boundary

             vol (Y)
        ________________
        |S_1,S_1,S_1,S_1|
        |               |
        |               |
S (X)   |               |
        |               |
        |               |
        |               |
        |               |
        |S_2,S_2,S_2,S_2|
        |_______________|

// ============================================================================
*/

class heston_equation
{

  private:
    heat_data_transform_2d_ptr heat_data_trans_cfg_;
    pde_discretization_config_2d_ptr discretization_cfg_;
    heston_boundary_transform_ptr heston_boundary_;
    splitting_method_config_ptr splitting_method_cfg_;
    grid_transform_config_2d_ptr grid_trans_cfg_; // this may be removed as it is not used later
    heat_implicit_solver_config_ptr solver_cfg_;
    std::map<std::string, double> solver_config_details_;

    explicit heston_equation() = delete;

    void initialize(heat_data_config_2d_ptr const &heat_data_cfg, grid_config_hints_2d_ptr const &grid_config_hints,
                    boundary_2d_ptr const &vertical_upper_boundary_ptr,
                    boundary_2d_pair const &horizontal_boundary_pair);

  public:
    explicit heston_equation(
        heat_data_config_2d_ptr const &heat_data_config, pde_discretization_config_2d_ptr const &discretization_config,
        boundary_2d_ptr const &vertical_upper_boundary_ptr, boundary_2d_pair const &horizontal_boundary_pair,
        splitting_method_config_ptr const &splitting_method_config, grid_config_hints_2d_ptr const &grid_config_hints,
        heat_implicit_solver_config_ptr const &solver_config =
            default_heat_solver_configs::host_fwd_dssolver_euler_solver_config_ptr,
        std::map<std::string, double> const &solver_config_details = std::map<std::string, double>());

    ~heston_equation();

    heston_equation(heston_equation const &) = delete;
    heston_equation(heston_equation &&) = delete;
    heston_equation &operator=(heston_equation const &) = delete;
    heston_equation &operator=(heston_equation &&) = delete;

    /**
     * Get the final solution of the PDE
     *
     * \param solution -  2D container for solution
     */
    LSS_API void solve(container_2d<by_enum::Row> &solution);

    /**
     * Get all solutions in time (surface) of the PDE
     *
     * \param solutions - 3D container for all the solutions in time
     */
    void solve(container_3d<by_enum::Row> &solutions);
};

} // namespace implicit_solvers

namespace explicit_solvers
{

/*!
============================================================================
Represents general variable coefficient Heston type equation

u_t = a(t,x,y)*u_xx + b(t,x,y)*u_yy + c(t,x,y)*u_xy + d(t,x,y)*u_x +
e(t,x,y)*u_y + f(t,x,y)*u + F(t,x,y)

t > 0, x_1 < x < x_2, y_1 < y < y_2

with initial condition:

u(0,x,y) = G(x,y)

or terminal condition:

u(T,x,y) = G(x,y)

horizontal_boundary_pair = S = (S_1,S_2) boundary

             vol (Y)
        ________________
        |S_1,S_1,S_1,S_1|
        |               |
        |               |
S (X)   |               |
        |               |
        |               |
        |               |
        |               |
        |S_2,S_2,S_2,S_2|
        |_______________|

// ============================================================================
*/

class heston_equation
{
  private:
    heat_data_transform_2d_ptr heat_data_trans_cfg_;
    pde_discretization_config_2d_ptr discretization_cfg_;
    heston_boundary_transform_ptr heston_boundary_;
    grid_transform_config_2d_ptr grid_trans_cfg_; // this may be removed as it is not used later
    heat_explicit_solver_config_ptr solver_cfg_;

    explicit heston_equation() = delete;

    void initialize(heat_data_config_2d_ptr const &heat_data_cfg, grid_config_hints_2d_ptr const &grid_config_hints,
                    boundary_2d_ptr const &vertical_upper_boundary_ptr,
                    boundary_2d_pair const &horizontal_boundary_pair);

  public:
    explicit heston_equation(heat_data_config_2d_ptr const &heat_data_config,
                             pde_discretization_config_2d_ptr const &discretization_config,
                             boundary_2d_ptr const &vertical_upper_boundary_ptr,
                             boundary_2d_pair const &horizontal_boundary_pair,
                             grid_config_hints_2d_ptr const &grid_config_hints,
                             heat_explicit_solver_config_ptr const &solver_config =
                                 default_heat_solver_configs::dev_expl_fwd_euler_solver_config_ptr);

    ~heston_equation();

    heston_equation(heston_equation const &) = delete;
    heston_equation(heston_equation &&) = delete;
    heston_equation &operator=(heston_equation const &) = delete;
    heston_equation &operator=(heston_equation &&) = delete;

    /**
     * Get the final solution of the PDE
     *
     * \param solution - 2D container for solution
     */
    LSS_API void solve(container_2d<by_enum::Row> &solution);

    /**
     * Get all solutions in time (surface) of the PDE
     *
     * \param solutions - 3D container for all the solutions in time
     */
    void solve(container_3d<by_enum::Row> &solutions);
};

} // namespace explicit_solvers

} // namespace two_dimensional

} // namespace lss_pde_solvers

#endif ///_LSS_HESTON_EQUATION_HPP_
