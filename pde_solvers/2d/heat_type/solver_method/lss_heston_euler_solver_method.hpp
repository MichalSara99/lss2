#if !defined(_LSS_HESTON_EULER_SOLVER_METHOD_HPP_)
#define _LSS_HESTON_EULER_SOLVER_METHOD_HPP_

#include "../../../../common/lss_enumerations.hpp"
#include "../../../../containers/lss_container_2d.hpp"
#include "../../../../discretization/lss_grid_config.hpp"
#include "../../../lss_pde_discretization_config.hpp"
#include "../explicit_coefficients/lss_heston_euler_coefficients.hpp"
#include "lss_heston_explicit_solver_method.hpp"

namespace lss_pde_solvers
{

using lss_containers::container_2d;
using lss_grids::grid_config_2d_ptr;

namespace two_dimensional
{

class explicit_heston_scheme
{

  public:
    static void rhs(heston_euler_coefficients_ptr const &cfs, grid_config_2d_ptr const &grid_config,
                    container_2d<by_enum::Row> const &input, double const &time, container_2d<by_enum::Row> &solution);

    static void rhs_source(heston_euler_coefficients_ptr const &cfs, grid_config_2d_ptr const &grid_config,
                           container_2d<by_enum::Row> const &input, double const &time,
                           container_2d<by_enum::Row> const &inhom_input, container_2d<by_enum::Row> &solution);
};

/**
    heston_euler_solver_method object
*/

class heston_euler_solver_method : public heston_explicit_solver_method
{

  private:
    // scheme coefficients:
    heston_euler_coefficients_ptr coefficients_;
    sptr_t<container_2d<by_enum::Row>> source_;

    explicit heston_euler_solver_method() = delete;

    void initialize(bool is_heat_source_set);

  public:
    explicit heston_euler_solver_method(heston_euler_coefficients_ptr const &coefficients,
                                        grid_config_2d_ptr const &grid_config, bool is_heat_source_set);

    ~heston_euler_solver_method();

    heston_euler_solver_method(heston_euler_solver_method const &) = delete;
    heston_euler_solver_method(heston_euler_solver_method &&) = delete;
    heston_euler_solver_method &operator=(heston_euler_solver_method const &) = delete;
    heston_euler_solver_method &operator=(heston_euler_solver_method &&) = delete;

    void solve(container_2d<by_enum::Row> &prev_solution, double const &time,
               container_2d<by_enum::Row> &solution) override;

    void solve(container_2d<by_enum::Row> &prev_solution, double const &time,
               std::function<double(double, double, double)> const &heat_source,
               container_2d<by_enum::Row> &solution) override;
};

} // namespace two_dimensional

} // namespace lss_pde_solvers

#endif ///_LSS_HESTON_EULER_SOLVER_METHOD_HPP_
