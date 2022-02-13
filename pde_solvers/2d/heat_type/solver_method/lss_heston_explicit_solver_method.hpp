#if !defined(_LSS_HESTON_EXPLICIT_SOLVER_METHOD_HPP_)
#define _LSS_HESTON_EXPLICIT_SOLVER_METHOD_HPP_

#include <functional>

#include "../../../../common/lss_enumerations.hpp"
#include "../../../../common/lss_utility.hpp"
#include "../../../../containers/lss_container_2d.hpp"
#include "../../../../discretization/lss_grid_config.hpp"

namespace lss_pde_solvers
{

using lss_containers::container_2d;
using lss_enumerations::by_enum;
using lss_grids::grid_config_2d_ptr;
using lss_utility::sptr_t;

namespace two_dimensional
{

class heston_explicit_solver_method
{

  protected:
    grid_config_2d_ptr grid_cfg_;

    explicit heston_explicit_solver_method() = delete;

  public:
    explicit heston_explicit_solver_method(grid_config_2d_ptr const &grid_config);

    virtual ~heston_explicit_solver_method();

    virtual void solve(container_2d<by_enum::Row> &prev_solution, double const &time,
                       container_2d<by_enum::Row> &solution) = 0;

    virtual void solve(container_2d<by_enum::Row> &prev_solution, double const &time,
                       std::function<double(double, double, double)> const &heat_source,
                       container_2d<by_enum::Row> &solution) = 0;
};

using heston_explicit_solver_method_ptr = sptr_t<heston_explicit_solver_method>;

} // namespace two_dimensional

} // namespace lss_pde_solvers

#endif ///_LSS_HESTON_EXPLICIT_SOLVER_METHOD_HPP_
