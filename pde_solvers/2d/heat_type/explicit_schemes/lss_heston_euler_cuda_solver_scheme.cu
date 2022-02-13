#include "lss_heston_euler_cuda_solver_scheme.hpp"
#include "../solver_method/lss_heston_euler_cuda_solver_method.hpp"

#include "../../../../common/lss_macros.hpp"

#define THREADS_PER_BLOCK_X 16
#define THREADS_PER_BLOCK_Y 16

namespace lss_pde_solvers
{

namespace two_dimensional
{

void heston_euler_cuda_kernel::launch(double time, thrust::device_vector<double> const &input,
                                      thrust::device_vector<double> &solution)
{
    discretize_coefficients(time);
    const unsigned int threads_per_block_x = THREADS_PER_BLOCK_X;
    const unsigned int threads_per_block_y = THREADS_PER_BLOCK_Y;
    const unsigned int blocks_per_grid_x =
        static_cast<unsigned int>(size_x_ + threads_per_block_x - 1) / threads_per_block_x;
    const unsigned int blocks_per_grid_y =
        static_cast<unsigned int>(size_y_ + threads_per_block_y - 1) / threads_per_block_y;
    const dim3 blocks_per_grid(blocks_per_grid_y, blocks_per_grid_x);
    const dim3 threads_per_block(threads_per_block_y, threads_per_block_x);
    double *raw_m = thrust::raw_pointer_cast(d_m_.data());
    double *raw_m_tilde = thrust::raw_pointer_cast(d_m_tilde_.data());
    double *raw_p = thrust::raw_pointer_cast(d_p_.data());
    double *raw_p_tilde = thrust::raw_pointer_cast(d_p_tilde_.data());
    double *raw_c = thrust::raw_pointer_cast(d_c_.data());
    double *raw_z = thrust::raw_pointer_cast(d_z_.data());
    double *raw_w = thrust::raw_pointer_cast(d_w_.data());
    const double *raw_input = thrust::raw_pointer_cast(input.data());
    double *raw_solution = thrust::raw_pointer_cast(solution.data());
    heston_core_kernel<<<blocks_per_grid, threads_per_block>>>(raw_m, raw_m_tilde, raw_p, raw_p_tilde, raw_c, raw_z,
                                                               raw_w, raw_input, raw_solution, size_x_, size_y_);
}

void heston_euler_cuda_kernel::launch(double time, thrust::device_vector<double> const &input,
                                      thrust::device_vector<double> const &source,
                                      thrust::device_vector<double> &solution)
{
    discretize_coefficients(time);
    const unsigned int threads_per_block_x = THREADS_PER_BLOCK_X;
    const unsigned int threads_per_block_y = THREADS_PER_BLOCK_Y;
    const unsigned int blocks_per_grid_x =
        static_cast<unsigned int>(size_x_ + threads_per_block_x - 1) / threads_per_block_x;
    const unsigned int blocks_per_grid_y =
        static_cast<unsigned int>(size_y_ + threads_per_block_y - 1) / threads_per_block_y;
    const dim3 blocks_per_grid(blocks_per_grid_y, blocks_per_grid_x);
    const dim3 threads_per_block(threads_per_block_y, threads_per_block_x);
    double *raw_m = thrust::raw_pointer_cast(d_m_.data());
    double *raw_m_tilde = thrust::raw_pointer_cast(d_m_tilde_.data());
    double *raw_p = thrust::raw_pointer_cast(d_p_.data());
    double *raw_p_tilde = thrust::raw_pointer_cast(d_p_tilde_.data());
    double *raw_c = thrust::raw_pointer_cast(d_c_.data());
    double *raw_z = thrust::raw_pointer_cast(d_z_.data());
    double *raw_w = thrust::raw_pointer_cast(d_w_.data());
    const double *raw_source = thrust::raw_pointer_cast(source.data());
    const double *raw_input = thrust::raw_pointer_cast(input.data());
    double *raw_solution = thrust::raw_pointer_cast(solution.data());
    heston_core_kernel<<<blocks_per_grid, threads_per_block>>>(raw_m, raw_m_tilde, raw_p, raw_p_tilde, raw_c, raw_z,
                                                               raw_w, k_, raw_input, raw_source, raw_solution, size_x_,
                                                               size_y_);
}

bool heston_euler_cuda_scheme::is_stable(heston_implicit_coefficients_ptr const &coefficients)
{
    // TODO: this needs to be implemented !!!
    return true;
}

void heston_euler_cuda_scheme::initialize(heston_implicit_coefficients_ptr const &coefficients)
{
    LSS_ASSERT(is_stable(coefficients) == true, "The chosen scheme is not stable");
    euler_coeffs_ = std::make_shared<heston_euler_coefficients>(coefficients);
    heston_boundary_ = std::make_shared<heston_explicit_boundary_solver>(coefficients, grid_cfg_);
}

heston_euler_cuda_scheme::heston_euler_cuda_scheme(heston_implicit_coefficients_ptr const &coefficients,
                                                   boundary_2d_ptr const &vertical_upper_boundary_ptr,
                                                   boundary_2d_pair const &horizontal_boundary_pair,
                                                   pde_discretization_config_2d_ptr const &discretization_config,
                                                   grid_config_2d_ptr const &grid_config)
    : boundary_ver_{vertical_upper_boundary_ptr}, boundary_pair_hor_{horizontal_boundary_pair},
      discretization_cfg_{discretization_config}, grid_cfg_{grid_config}
{
    initialize(coefficients);
}

heston_euler_cuda_scheme::~heston_euler_cuda_scheme()
{
}

} // namespace two_dimensional

} // namespace lss_pde_solvers
