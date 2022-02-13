#include "lss_heston_euler_cuda_solver_method.hpp"

#include "../../../../discretization/lss_grid.hpp"
#include "../../../../discretization/lss_grid_config.hpp"

namespace lss_pde_solvers
{
namespace two_dimensional
{

__global__ void heston_core_kernel(double const *m_coeff, double const *m_tilde_coeff, double const *p_coeff,
                                   double const *p_tilde_coeff, double const *c_coeff, double const *z_coeff,
                                   double const *w_coeff, double const *input, double *solution,
                                   const std::size_t size_x, const std::size_t size_y)
{
    const std::size_t c_id = blockDim.x * blockIdx.x + threadIdx.x; // size_y
    const std::size_t r_id = blockDim.y * blockIdx.y + threadIdx.y; // size_x
    const std::size_t tid = c_id + r_id * size_y;
    if (c_id == 0)
        return;
    if (c_id >= (size_y - 1))
        return;
    if (r_id == 0)
        return;
    if (r_id >= (size_x - 1))
        return;
    // cross neighbours:
    const std::size_t up_tid = tid - size_y;
    const std::size_t down_tid = tid + size_y;
    const std::size_t right_tid = tid + 1;
    const std::size_t left_tid = tid - 1;
    // star neighbours:
    const std::size_t up_r_tid = tid - size_y + 1;
    const std::size_t up_l_tid = tid - size_y - 1;
    const std::size_t down_r_tid = tid + size_y + 1;
    const std::size_t down_l_tid = tid + size_y - 1;
    const double one = 1.0;

    solution[tid] = (c_coeff[tid] * input[up_l_tid]) + (m_coeff[tid] * input[up_tid]) -
                    (c_coeff[tid] * input[up_r_tid]) + (m_tilde_coeff[tid] * input[left_tid]) +
                    ((one - z_coeff[tid] - w_coeff[tid]) * input[tid]) + (p_tilde_coeff[tid] * input[right_tid]) -
                    (c_coeff[tid] * input[down_l_tid]) + (p_coeff[tid] * input[down_tid]) +
                    (c_coeff[tid] * input[down_r_tid]);
}

__global__ void heston_core_kernel(double const *m_coeff, double const *m_tilde_coeff, double const *p_coeff,
                                   double const *p_tilde_coeff, double const *c_coeff, double const *z_coeff,
                                   double const *w_coeff, double time_step, double const *input, double const *source,
                                   double *solution, const std::size_t size_x, const std::size_t size_y)
{
    const std::size_t c_id = blockDim.x * blockIdx.x + threadIdx.x; // size_y
    const std::size_t r_id = blockDim.y * blockIdx.y + threadIdx.y; // size_x
    const std::size_t tid = c_id + r_id * size_y;
    if (c_id == 0)
        return;
    if (c_id >= (size_y - 1))
        return;
    if (r_id == 0)
        return;
    if (r_id >= (size_x - 1))
        return;
    // cross neighbours:
    const std::size_t up_tid = tid - size_y;
    const std::size_t down_tid = tid + size_y;
    const std::size_t right_tid = tid + 1;
    const std::size_t left_tid = tid - 1;
    // star neighbours:
    const std::size_t up_r_tid = tid - size_y + 1;
    const std::size_t up_l_tid = tid - size_y - 1;
    const std::size_t down_r_tid = tid + size_y + 1;
    const std::size_t down_l_tid = tid + size_y - 1;
    const double one = 1.0;

    solution[tid] = (c_coeff[tid] * input[up_l_tid]) + (m_coeff[tid] * input[up_tid]) -
                    (c_coeff[tid] * input[up_r_tid]) + (m_tilde_coeff[tid] * input[left_tid]) +
                    ((one - z_coeff[tid] - w_coeff[tid]) * input[tid]) + (p_tilde_coeff[tid] * input[right_tid]) -
                    (c_coeff[tid] * input[down_l_tid]) + (p_coeff[tid] * input[down_tid]) +
                    (c_coeff[tid] * input[down_r_tid]) + (time_step * source[tid]);
}

void heston_euler_cuda_kernel::initialize(heston_euler_coefficients_ptr const &coefficients)
{
    size_x_ = coefficients->space_size_x_;
    size_y_ = coefficients->space_size_y_;
    const std::size_t total_size = size_x_ * size_y_;
    k_ = coefficients->k_;
    m_ = coefficients->M_;
    m_tilde_ = coefficients->M_tilde_;
    p_ = coefficients->P_;
    p_tilde_ = coefficients->P_tilde_;
    c_ = coefficients->C_;
    z_ = coefficients->Z_;
    w_ = coefficients->W_;
    // on host:
    h_m_.resize(total_size);
    h_m_tilde_.resize(total_size);
    h_p_.resize(total_size);
    h_p_tilde_.resize(total_size);
    h_c_.resize(total_size);
    h_z_.resize(total_size);
    h_w_.resize(total_size);
    // on device:
    d_m_.resize(total_size);
    d_m_tilde_.resize(total_size);
    d_p_.resize(total_size);
    d_p_tilde_.resize(total_size);
    d_c_.resize(total_size);
    d_z_.resize(total_size);
    d_w_.resize(total_size);
}

void heston_euler_cuda_kernel::discretize_coefficients(double time)
{
    // discretize on host
    d_2d::of_function(grid_cfg_, time, m_, size_x_, size_y_, h_m_);
    d_2d::of_function(grid_cfg_, time, m_tilde_, size_x_, size_y_, h_m_tilde_);
    d_2d::of_function(grid_cfg_, time, p_, size_x_, size_y_, h_p_);
    d_2d::of_function(grid_cfg_, time, p_tilde_, size_x_, size_y_, h_p_tilde_);
    d_2d::of_function(grid_cfg_, time, c_, size_x_, size_y_, h_c_);
    d_2d::of_function(grid_cfg_, time, z_, size_x_, size_y_, h_z_);
    d_2d::of_function(grid_cfg_, time, w_, size_x_, size_y_, h_w_);
    // copy to device
    thrust::copy(h_m_.begin(), h_m_.end(), d_m_.begin());
    thrust::copy(h_m_tilde_.begin(), h_m_tilde_.end(), d_m_tilde_.begin());
    thrust::copy(h_p_.begin(), h_p_.end(), d_p_.begin());
    thrust::copy(h_p_tilde_.begin(), h_p_tilde_.end(), d_p_tilde_.begin());
    thrust::copy(h_c_.begin(), h_c_.end(), d_c_.begin());
    thrust::copy(h_z_.begin(), h_z_.end(), d_z_.begin());
    thrust::copy(h_w_.begin(), h_w_.end(), d_w_.begin());
}

heston_euler_cuda_kernel::heston_euler_cuda_kernel(heston_euler_coefficients_ptr const &coefficients,
                                                   grid_config_2d_ptr const &grid_config)
    : grid_cfg_{grid_config}
{
    initialize(coefficients);
}

heston_euler_cuda_kernel::~heston_euler_cuda_kernel()
{
}

void explicit_heston_cuda_scheme::rhs(heston_euler_coefficients_ptr const &cfs,
                                      heston_euler_cuda_kernel_ptr const &kernel, grid_config_2d_ptr const &grid_config,
                                      container_2d<by_enum::Row> const &input, double const &time,
                                      container_2d<by_enum::Row> &solution)
{
    // light-weight object with cuda kernel computing the solution:
    thrust::device_vector<double> d_input(input.data());
    thrust::device_vector<double> d_solution(solution.data());
    kernel->launch(time, d_input, d_solution);
    container_t tmp_solution(solution.data().size());
    thrust::copy(d_solution.begin(), d_solution.end(), tmp_solution.begin());
    solution.from_data(tmp_solution);
}

void explicit_heston_cuda_scheme::rhs_source(heston_euler_coefficients_ptr const &cfs,
                                             heston_euler_cuda_kernel_ptr const &kernel,
                                             grid_config_2d_ptr const &grid_config,
                                             container_2d<by_enum::Row> const &input, double const &time,
                                             container_2d<by_enum::Row> const &inhom_input,
                                             container_2d<by_enum::Row> &solution)
{
    // light-weight object with cuda kernel computing the solution:
    thrust::device_vector<double> d_input(input.data());
    thrust::device_vector<double> d_inhom_input(inhom_input.data());
    thrust::device_vector<double> d_solution(solution.data());
    kernel->launch(time, d_input, d_inhom_input, d_solution);
    container_t tmp_solution(solution.data().size());
    thrust::copy(d_solution.begin(), d_solution.end(), tmp_solution.begin());
    solution.from_data(tmp_solution);
}

void heston_euler_cuda_solver_method::initialize(bool is_heat_source_set)
{
    kernel_ = std::make_shared<heston_euler_cuda_kernel>(coefficients_, grid_cfg_);
    if (is_heat_source_set)
    {
        source_ =
            std::make_shared<container_2d<by_enum::Row>>(coefficients_->space_size_x_, coefficients_->space_size_y_);
    }
}

heston_euler_cuda_solver_method::heston_euler_cuda_solver_method(heston_euler_coefficients_ptr const &coefficients,
                                                                 grid_config_2d_ptr const &grid_config,
                                                                 bool is_heat_source_set)
    : coefficients_{coefficients}, heston_explicit_solver_method(grid_config){

    initialize(is_heat_source_set);
}

heston_euler_cuda_solver_method::~heston_euler_cuda_solver_method()
{
}

void heston_euler_cuda_solver_method::solve(container_2d<by_enum::Row> &prev_solution, double const &time,
                                            container_2d<by_enum::Row> &solution)
{
    explicit_heston_cuda_scheme::rhs(coefficients_, kernel_, grid_cfg_, prev_solution, time, solution);
}

void heston_euler_cuda_solver_method::solve(container_2d<by_enum::Row> &prev_solution, double const &time,
                                            std::function<double(double, double, double)> const &heat_source,
                                            container_2d<by_enum::Row> &solution)
{
    d_2d::of_function(grid_cfg_, time, heat_source, *source_);
    explicit_heston_cuda_scheme::rhs_source(coefficients_, kernel_, grid_cfg_, prev_solution, time, *source_, solution);
}

} // namespace two_dimensional

} // namespace lss_pde_solvers
