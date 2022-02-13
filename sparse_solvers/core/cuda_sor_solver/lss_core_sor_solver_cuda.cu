#include "lss_core_sor_solver_cuda.hpp"

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

namespace lss_core_sor_solver
{

 core_sor_solver_cuda::core_sor_solver_cuda(){};


core_sor_solver_cuda::core_sor_solver_cuda(std::size_t system_size)
    : system_size_{system_size}
{
}

 core_sor_solver_cuda::~core_sor_solver_cuda()
{
}

void core_sor_solver_cuda::set_omega(double value)
{
   LSS_ASSERT(
      (value > static_cast<double>(0.0)) && (value < static_cast<double>(2.0)),
               "relaxation parameter must be inside (0,2) range");
    omega_ = value;
}

void core_sor_solver_cuda::set_rhs(container_t const &rhs)
{
    LSS_ASSERT(rhs.size() == system_size_, "Inncorect size for right-hand side");
    b_ = rhs;
}


void core_sor_solver_cuda::set_flat_sparse_matrix(flat_matrix matrix)
{
    LSS_ASSERT(matrix.rows() == system_size_, "Inncorect number of rows for the flat_raw_matrix");
    LSS_ASSERT(matrix.columns() == system_size_, "Inncorect number of columns for the flat_raw_matrix");
    matrix_data_ptr_ = std::make_unique<flat_matrix>(std::move(matrix));
}


 bool core_sor_solver_cuda::is_diagonally_dominant()
{
    LSS_VERIFY(matrix_data_ptr_, "flat_matrix has not been provided.");
    LSS_ASSERT(b_.size() == system_size_, "Incorrect size for right-hand side");
    // first make sure the matrix is row-major sorted
    matrix_data_ptr_->sort(flat_matrix_sort_enum::RowMajor);
    double diag{};
    std::tuple<std::size_t, std::size_t, double> non_diag{};
    double sum{};
    std::size_t cols{};
    // for index of the flat_matrix element:
    std::size_t flt_idx{0};
    for (std::size_t r = 0; r < matrix_data_ptr_->rows(); ++r, flt_idx += cols)
    {
      sum = static_cast<double>(0.0);
        diag = std::abs(matrix_data_ptr_->diagonal_at_row(r));
        cols = matrix_data_ptr_->non_zero_column_size(r);
        for (std::size_t c = flt_idx; c < cols + flt_idx; ++c)
        {
            non_diag = matrix_data_ptr_->at(c);
            if (std::get<0>(non_diag) != std::get<1>(non_diag))
                sum += std::abs(std::get<2>(non_diag));
        }
        if (diag < sum)
            return false;
    }
    return true;
}


void core_sor_solver_cuda::kernel(container_t &solution)
{
    LSS_ASSERT(is_diagonally_dominant() == true, "flat_raw_matrix isd not diagonally dominant.");
    // set initial step:
    std::size_t step{0};
    // set iter_limit:
    std::size_t iter_limit = sor_solver_cuda_traits<double>::iteration_limit();
    // set tolerance:
    double const tol = sor_solver_cuda_traits<double>::tolerance();
    // for error:
    double error{};
    // make copies to pass to sor_solver_core:
    std::size_t total_size{matrix_data_ptr_->size()};

    thrust::host_vector<double> h(b_.size());
    copy(h, b_);
    thrust::device_vector<double> rhs_values(h);
    thrust::device_vector<double> matrix_values(total_size);
    thrust::device_vector<std::size_t> non_zero_column_idx(total_size);
    std::tuple<std::size_t, std::size_t, double> triplet{};
    for (std::size_t t = 0; t < total_size; ++t)
    {
        triplet = matrix_data_ptr_->at(t);
        non_zero_column_idx[t] = std::get<1>(triplet);
        matrix_values[t] = std::get<2>(triplet);
    }
    std::size_t row_size{matrix_data_ptr_->rows()};
    thrust::device_vector<double> diagonal_values(row_size);
    thrust::device_vector<std::size_t> row_start_idx(row_size);
    thrust::device_vector<std::size_t> row_end_idx(row_size);
    std::size_t end_cnt{0};
    std::size_t non_zero{0};
    for (std::size_t r = 0; r < row_size; ++r)
    {
        diagonal_values[r] = matrix_data_ptr_->diagonal_at_row(r);
        row_start_idx[r] = end_cnt;
        non_zero = matrix_data_ptr_->non_zero_column_size(r);
        end_cnt += non_zero;
        row_end_idx[r] = end_cnt - 1;
    }
    // initialize sor_solver_core object:
    sor_solver_core core_solver(matrix_values, non_zero_column_idx, row_size, rhs_values, diagonal_values,
                                       row_start_idx, row_end_idx);

    thrust::host_vector<double> hsol(solution.size());
    copy(hsol, solution);
    thrust::device_vector<double> sol(hsol);
    thrust::device_vector<double> new_sol(solution.size());
    thrust::device_vector<double> errors(solution.size());

    while (iter_limit > step)
    {
        core_solver(sol, new_sol, errors, omega_);
        error = thrust::reduce(errors.begin(), errors.end(),
                             decltype(errors)::value_type(0.0),
                             thrust::plus<double>());
        if (error <= tol)
            break;
        sol = new_sol;
        step++;
    }
    copy(solution,sol);
}


void core_sor_solver_cuda::solve(container_t &solution)
{
    LSS_ASSERT(solution.size() == system_size_, "Incorrect size of solution container");
    kernel(solution);
}

container_t const core_sor_solver_cuda::solve()
{
    container_t solution(system_size_);
    kernel(solution);
    return solution;
}

// double precision specialization:



} // namespace lss_core_sor_solver
