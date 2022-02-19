#include "lss_heat_euler_cuda_scheme_2d.hpp"
#include "../solver_method/lss_heat_euler_cuda_solver_method_2d.hpp"

#include "../../../../common/lss_macros.hpp"

#define THREADS_PER_BLOCK_X 16
#define THREADS_PER_BLOCK_Y 16

namespace lss_pde_solvers
{

namespace two_dimensional
{

void heat_euler_cuda_kernel_2d::launch(double time, thrust::device_vector<double> const &input,
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
    heat_core_kernel_2d<<<blocks_per_grid, threads_per_block>>>(raw_m, raw_m_tilde, raw_p, raw_p_tilde, raw_c, raw_z,
                                                               raw_w, raw_input, raw_solution, size_x_, size_y_);
}

void heat_euler_cuda_kernel_2d::launch(double time, thrust::device_vector<double> const &input,
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
    heat_core_kernel_2d<<<blocks_per_grid, threads_per_block>>>(raw_m, raw_m_tilde, raw_p, raw_p_tilde, raw_c, raw_z,
                                                               raw_w, k_, raw_input, raw_source, raw_solution, size_x_,
                                                               size_y_);
}

bool heat_euler_cuda_scheme_2d::is_stable(heat_coefficients_2d_ptr const &coefficients) {
    // TODO: this needs to be implemented !!!
    return true;
}

void heat_euler_cuda_scheme_2d::initialize(heat_coefficients_2d_ptr const &coefficients) {
    LSS_ASSERT(is_stable(coefficients) == true, "The chosen scheme is not stable");
    euler_coeffs_ = std::make_shared<heat_euler_coefficients_2d>(coefficients);
    heston_boundary_ = std::make_shared<heston_boundary_solver>(coefficients, grid_cfg_);
}

heat_euler_cuda_scheme_2d::heat_euler_cuda_scheme_2d(heat_coefficients_2d_ptr const &coefficients,
                                                   boundary_2d_ptr const &vertical_upper_boundary_ptr,
                                                   boundary_2d_pair const &horizontal_boundary_pair,
                                                   pde_discretization_config_2d_ptr const &discretization_config,
                                                   grid_config_2d_ptr const &grid_config)
    : boundary_ver_{vertical_upper_boundary_ptr}, boundary_pair_hor_{horizontal_boundary_pair},
      discretization_cfg_{discretization_config}, grid_cfg_{grid_config}
{
    initialize(coefficients);
}

heat_euler_cuda_scheme_2d::~heat_euler_cuda_scheme_2d() {}

void heat_euler_cuda_scheme_2d::operator()(matrix_2d &prev_solution, matrix_2d &next_solution, bool is_heat_sourse_set,
                                          std::function<double(double, double, double)> const &heat_source,
                                          traverse_direction_enum traverse_dir) {
    auto const &timer = discretization_cfg_->time_range();
    const double k = discretization_cfg_->time_step();
    // last time index:
    const std::size_t last_time_idx = discretization_cfg_->number_of_time_points() - 1;
    auto const &solver_method_ptr =
        std::make_shared<heat_euler_cuda_solver_method_2d>(euler_coeffs_, grid_cfg_, is_heat_sourse_set);
    if (is_heat_sourse_set) {
      throw std::exception("Not yet implemented");
        // TODO: to be implemented!!!
        //  loop::run(solver_method_ptr, boundary_pair_, timer, last_time_idx, k,
        //  traverse_dir, heat_source, solution);
    } else {
        explicit_time_loop_2d::run(solver_method_ptr, heston_boundary_,
                                   boundary_pair_hor_, boundary_ver_, grid_cfg_,
                                   timer, last_time_idx, k, traverse_dir,
                                   prev_solution, next_solution);
    }
}

void heat_euler_cuda_scheme_2d::operator()(matrix_2d &prev_solution, matrix_2d &next_solution, bool is_heat_sourse_set,
                                            std::function<double(double, double, double)> const &heat_source,
                                            traverse_direction_enum traverse_dir, matrix_3d &solutions) {
  auto const &timer = discretization_cfg_->time_range();
  const double k = discretization_cfg_->time_step();
  // last time index:
  const std::size_t last_time_idx =
      discretization_cfg_->number_of_time_points() - 1;
  auto const &solver_method_ptr =
      std::make_shared<heat_euler_cuda_solver_method_2d>(
          euler_coeffs_, grid_cfg_, is_heat_sourse_set);
  if (is_heat_sourse_set) {
        throw std::exception("Not yet implemented");
        // TODO: to be implemented!!!
        //  loop::run(solver_method_ptr, boundary_pair_, timer, last_time_idx, k,
        //  traverse_dir, heat_source, solution);
  } else {
    explicit_time_loop_2d::run_with_stepping(solver_method_ptr, heston_boundary_,
                               boundary_pair_hor_, boundary_ver_, grid_cfg_,
                               timer, last_time_idx, k, traverse_dir,
                               prev_solution, next_solution, solutions);
  }
}

} // namespace two_dimensional

} // namespace lss_pde_solvers
