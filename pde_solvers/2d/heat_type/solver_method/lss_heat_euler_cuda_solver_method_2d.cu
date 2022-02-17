#include "lss_heat_euler_cuda_solver_method_2d.hpp"

#include "../../../../discretization/lss_grid.hpp"
#include "../../../../discretization/lss_grid_config.hpp"
#include"../../../../common/lss_utility.hpp"

namespace lss_pde_solvers
{
namespace two_dimensional
{

    using lss_utility::copy;

__global__ void heat_core_kernel_2d(double const *m_coeff, double const *m_tilde_coeff, double const *p_coeff,
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

__global__ void heat_core_kernel_2d(double const *m_coeff, double const *m_tilde_coeff, double const *p_coeff,
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

void heat_euler_cuda_kernel_2d::initialize(heat_euler_coefficients_2d_ptr const &coefficients) {
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

void heat_euler_cuda_kernel_2d::discretize_coefficients(double time) {
    // discretize on host
    d_2d::of_function(grid_cfg_, time, m_, size_x_, size_y_, h_m_);
    d_2d::of_function(grid_cfg_, time, m_tilde_, size_x_, size_y_, h_m_tilde_);
    d_2d::of_function(grid_cfg_, time, p_, size_x_, size_y_, h_p_);
    d_2d::of_function(grid_cfg_, time, p_tilde_, size_x_, size_y_, h_p_tilde_);
    d_2d::of_function(grid_cfg_, time, c_, size_x_, size_y_, h_c_);
    d_2d::of_function(grid_cfg_, time, z_, size_x_, size_y_, h_z_);
    d_2d::of_function(grid_cfg_, time, w_, size_x_, size_y_, h_w_);
    // copy to host vector first:
    thrust::host_vector<double> m(h_m_.size());
    thrust::host_vector<double> m_tilde(h_m_tilde_.size());
    thrust::host_vector<double> p(h_p_.size());
    thrust::host_vector<double> p_tilde(h_p_tilde_.size());
    thrust::host_vector<double> c(h_c_.size());
    thrust::host_vector<double> z(h_z_.size());
    thrust::host_vector<double> w(h_w_.size());
    copy(m, h_m_);
    copy(m_tilde, h_m_tilde_);
    copy(p, h_p_);
    copy(p_tilde, h_p_tilde_);
    copy(c, h_c_);
    copy(z, h_z_);
    copy(w, h_w_);
    // copy to device
    thrust::copy(m.begin(), m.end(), d_m_.begin());
    thrust::copy(m_tilde.begin(), m_tilde.end(), d_m_tilde_.begin());
    thrust::copy(p.begin(), p.end(), d_p_.begin());
    thrust::copy(p_tilde.begin(), p_tilde.end(), d_p_tilde_.begin());
    thrust::copy(c.begin(), c.end(), d_c_.begin());
    thrust::copy(z.begin(), z.end(), d_z_.begin());
    thrust::copy(w.begin(), w.end(), d_w_.begin());
}

heat_euler_cuda_kernel_2d::heat_euler_cuda_kernel_2d(heat_euler_coefficients_2d_ptr const &coefficients,
                                                   grid_config_2d_ptr const &grid_config)
    : grid_cfg_{grid_config}
{
    initialize(coefficients);
}

heat_euler_cuda_kernel_2d::~heat_euler_cuda_kernel_2d() {}

void explicit_euler_cuda_rhs::rhs(heat_euler_coefficients_2d_ptr const &cfs,
                                  heat_euler_cuda_kernel_2d_ptr const &kernel,
                                  grid_config_2d_ptr const &grid_config,
                                  matrix_2d const &input, double const &time,
                                  matrix_2d &solution) {
    // light-weight object with cuda kernel computing the solution:
    thrust::host_vector<double> h_input(input.total_size());
    thrust::host_vector<double> h_solution(solution.total_size());
    copy(h_input, input.data());
    copy(h_solution, solution.data());
    thrust::device_vector<double> d_input(h_input);
    thrust::device_vector<double> d_solution(h_solution);
    kernel->launch(time, d_input, d_solution);
    thrust::copy(d_solution.begin(), d_solution.end(), h_solution.begin());
    container_t tmp_solution(solution.total_size());
    copy(tmp_solution, h_solution);
    solution.from_data(std::move(tmp_solution));
}

void explicit_euler_cuda_rhs::rhs_source(heat_euler_coefficients_2d_ptr const &cfs,
                                        heat_euler_cuda_kernel_2d_ptr const &kernel,
                                        grid_config_2d_ptr const &grid_config, 
                                        matrix_2d const &input, double const &time,
                                        matrix_2d const &inhom_input, matrix_2d &solution) {
    // light-weight object with cuda kernel computing the solution:
    thrust::host_vector<double> h_input(input.total_size());
    thrust::host_vector<double> h_inhom_input(inhom_input.total_size());
    thrust::host_vector<double> h_solution(solution.total_size());
    copy(h_input, input.data());
    copy(h_inhom_input, inhom_input.data());
    copy(h_solution, solution.data());
    thrust::device_vector<double> d_input(h_input);
    thrust::device_vector<double> d_inhom_input(h_inhom_input);
    thrust::device_vector<double> d_solution(h_solution);
    kernel->launch(time, d_input, d_inhom_input, d_solution);
    thrust::copy(d_solution.begin(), d_solution.end(), h_solution.begin());
    container_t tmp_solution(solution.total_size());
    copy(tmp_solution, h_solution);
    solution.from_data(std::move(tmp_solution));
}

void heat_euler_cuda_solver_method_2d::initialize(bool is_heat_source_set) {
    kernel_ = std::make_shared<heat_euler_cuda_kernel_2d>(coefficients_, grid_cfg_);
    if (is_heat_source_set)
    {
        source_ =
            std::make_shared<matrix_2d>(coefficients_->space_size_x_, coefficients_->space_size_y_);
    }
}

heat_euler_cuda_solver_method_2d::heat_euler_cuda_solver_method_2d(
    heat_euler_coefficients_2d_ptr const &coefficients,
                                                                        grid_config_2d_ptr const &grid_config,bool is_heat_source_set)
    : coefficients_{coefficients}, heat_explicit_solver_method_2d(grid_config) {

    initialize(is_heat_source_set);
}

heat_euler_cuda_solver_method_2d::~heat_euler_cuda_solver_method_2d() {}

void heat_euler_cuda_solver_method_2d::solve(matrix_2d &prev_solution,
                                             double const &time,
                                                matrix_2d &solution) {
  explicit_euler_cuda_rhs::rhs(coefficients_, kernel_, grid_cfg_, prev_solution,time, solution);
}

void heat_euler_cuda_solver_method_2d::solve(
    matrix_2d &prev_solution, double const &time,
                                                std::function<double(double, double, double)> const &heat_source,
                                                matrix_2d &solution) {
    d_2d::of_function(grid_cfg_, time, heat_source, *source_);
    explicit_euler_cuda_rhs::rhs_source(coefficients_, kernel_, grid_cfg_,prev_solution, time, *source_, solution);
}

} // namespace two_dimensional

} // namespace lss_pde_solvers
