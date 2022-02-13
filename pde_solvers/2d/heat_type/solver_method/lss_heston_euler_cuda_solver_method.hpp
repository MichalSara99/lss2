#if !defined(_LSS_HESTON_EULER_CUDA_SOLVER_METHOD_HPP_)
#define _LSS_HESTON_EULER_CUDA_SOLVER_METHOD_HPP_

#include <device_launch_parameters.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

#include "../../../../common/lss_enumerations.hpp"
#include "../../../../containers/lss_container_2d.hpp"
#include "../../../../discretization/lss_discretization.hpp"
#include "../../../../discretization/lss_grid_config.hpp"
#include "../explicit_coefficients/lss_heston_euler_coefficients.hpp"
#include "lss_heston_explicit_solver_method.hpp"

namespace lss_pde_solvers
{
namespace two_dimensional
{

using lss_containers::container_2d;
using lss_enumerations::by_enum;
using lss_grids::grid_config_2d_ptr;
using lss_utility::container_t;
using lss_utility::NaN;
using lss_utility::range;
using lss_utility::sptr_t;

using d_2d = discretization_2d<std::vector, std::allocator<double>>;

__global__ extern void heston_core_kernel(double const *m_coeff, double const *m_tilde_coeff, double const *p_coeff,
                                          double const *p_tilde_coeff, double const *c_coeff, double const *z_coeff,
                                          double const *w_coeff, double const *input, double *solution,
                                          const std::size_t size_x, const std::size_t size_y);

__global__ extern void heston_core_kernel(double const *m_coeff, double const *m_tilde_coeff, double const *p_coeff,
                                          double const *p_tilde_coeff, double const *c_coeff, double const *z_coeff,
                                          double const *w_coeff, double time_step, double const *input,
                                          double const *source, double *solution, const std::size_t size_x,
                                          const std::size_t size_y);

/**
 * heston_euler_cuda_kernel object
 */
class heston_euler_cuda_kernel
{

  private:
    double k_;
    std::size_t size_x_, size_y_;
    // device containers:
    thrust::device_vector<double> d_m_;
    thrust::device_vector<double> d_m_tilde_;
    thrust::device_vector<double> d_p_;
    thrust::device_vector<double> d_p_tilde_;
    thrust::device_vector<double> d_c_;
    thrust::device_vector<double> d_z_;
    thrust::device_vector<double> d_w_;
    // host containers:
    container_t h_m_;
    container_t h_m_tilde_;
    container_t h_p_;
    container_t h_p_tilde_;
    container_t h_c_;
    container_t h_z_;
    container_t h_w_;
    grid_config_2d_ptr grid_cfg_;
    // coefficients:
    std::function<double(double, double, double)> m_;
    std::function<double(double, double, double)> m_tilde_;
    std::function<double(double, double, double)> p_;
    std::function<double(double, double, double)> p_tilde_;
    std::function<double(double, double, double)> c_;
    std::function<double(double, double, double)> z_;
    std::function<double(double, double, double)> w_;

    void initialize(heston_euler_coefficients_ptr const &coefficients);

    void discretize_coefficients(double time);

  public:
    explicit heston_euler_cuda_kernel(heston_euler_coefficients_ptr const &coefficients,
                                      grid_config_2d_ptr const &grid_config);

    ~heston_euler_cuda_kernel();

    void launch(double time, thrust::device_vector<double> const &input, thrust::device_vector<double> &solution);

    void launch(double time, thrust::device_vector<double> const &input, thrust::device_vector<double> const &source,
                thrust::device_vector<double> &solution);
};

using heston_euler_cuda_kernel_ptr = sptr_t<heston_euler_cuda_kernel>;

class explicit_heston_cuda_scheme
{

  public:
    static void rhs(heston_euler_coefficients_ptr const &cfs, heston_euler_cuda_kernel_ptr const &kernel,
                    grid_config_2d_ptr const &grid_config, container_2d<by_enum::Row> const &input, double const &time,
                    container_2d<by_enum::Row> &solution);

    static void rhs_source(heston_euler_coefficients_ptr const &cfs, heston_euler_cuda_kernel_ptr const &kernel,
                           grid_config_2d_ptr const &grid_config, container_2d<by_enum::Row> const &input,
                           double const &time, container_2d<by_enum::Row> const &inhom_input,
                           container_2d<by_enum::Row> &solution);
};

/**
    heston_euler_cuda_solver_method object
*/

class heston_euler_cuda_solver_method : public heston_explicit_solver_method
{

  private:
    // scheme coefficients:
    heston_euler_coefficients_ptr coefficients_;
    // cuda kernel:
    heston_euler_cuda_kernel_ptr kernel_;
    sptr_t<container_2d<by_enum::Row>> source_;

    explicit heston_euler_cuda_solver_method() = delete;

    void initialize(bool is_heat_source_set);

  public:
    explicit heston_euler_cuda_solver_method(heston_euler_coefficients_ptr const &coefficients,
                                             grid_config_2d_ptr const &grid_config, bool is_heat_source_set);

    ~heston_euler_cuda_solver_method();

    heston_euler_cuda_solver_method(heston_euler_cuda_solver_method const &) = delete;
    heston_euler_cuda_solver_method(heston_euler_cuda_solver_method &&) = delete;
    heston_euler_cuda_solver_method &operator=(heston_euler_cuda_solver_method const &) = delete;
    heston_euler_cuda_solver_method &operator=(heston_euler_cuda_solver_method &&) = delete;

    void solve(container_2d<by_enum::Row> &prev_solution, double const &time,
               container_2d<by_enum::Row> &solution) override;

    void solve(container_2d<by_enum::Row> &prev_solution, double const &time,
               std::function<double(double, double, double)> const &heat_source,
               container_2d<by_enum::Row> &solution) override;
};

} // namespace two_dimensional

} // namespace lss_pde_solvers

#endif ///_LSS_HESTON_EULER_CUDA_SOLVER_METHOD_HPP_
