#include "lss_core_sor_solver.hpp"

#include "../traits/lss_sor_solver_traits.hpp"

namespace lss_core_sor_solver
{

using lss_containers::flat_matrix_sort_enum;
using lss_sor_solver_traits::sor_solver_traits;

core_sor_solver::core_sor_solver()
{
}

core_sor_solver::core_sor_solver(std::size_t system_size) : system_size_{system_size}
{
}

core_sor_solver::~core_sor_solver()
{
}

void core_sor_solver::set_omega(double value)
{
    LSS_ASSERT((value > static_cast<double>(0.0)) && (value < static_cast<double>(2.0)),
               "relaxation parameter must be inside (0,2) range");
    omega_ = value;
}

void lss_core_sor_solver::core_sor_solver::set_rhs(container_t const &rhs)
{
    LSS_ASSERT(rhs.size() == system_size_, "Inncorect size for right-hand side");
    b_ = rhs;
}

void core_sor_solver::set_flat_sparse_matrix(flat_matrix matrix)
{
    LSS_ASSERT(matrix.rows() == system_size_, "Inncorect number of rows for the flat_raw_matrix");
    LSS_ASSERT(matrix.columns() == system_size_, "Inncorect number of columns for the flat_raw_matrix");
    matrix_data_ptr_ = std::make_unique<flat_matrix>(std::move(matrix));
}

bool core_sor_solver::is_diagonally_dominant()
{
    LSS_VERIFY(matrix_data_ptr_, "flat_matrix has not been provided.");
    LSS_ASSERT(b_.size() == system_size_, "Incorrect size for right-hand side");
    // first make sure the matrix is row-major sorted
    matrix_data_ptr_->sort(flat_matrix_sort_enum::RowMajor);
    double diag{};
    std::tuple<std::size_t, std::size_t, double> non_diag;
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

void core_sor_solver::kernel(container_t &solution)
{
    LSS_ASSERT(is_diagonally_dominant() == true, "flat_raw_matrix isd not diagonally dominant.");
    // set initial step:
    std::size_t step{0};
    // set iter_limit:
    const std::size_t iter_limit = sor_solver_traits<double>::iteration_limit();
    // set tolerance:
    const double tol = sor_solver_traits<double>::tolerance();
    // for error:
    double error{};
    // for sigma_1,sigma_2:
    double sigma_1{};
    double sigma_2{};
    // for diagonal value:
    double diag{};
    // for new solution:
    container_t x_new(solution);
    // for number of columns:
    std::size_t cols{};
    // for flat_matrix element:
    std::tuple<std::size_t, std::size_t, double> mat_elm;
    // for flat_matrix row and column index of the element:
    std::size_t mat_r{};
    std::size_t mat_c{};
    // for flat_matrix element value:
    double mat_val{};
    // for index of the flat_matrix element:
    std::size_t flt_idx{0};
    const double one = static_cast<double>(1.0);

    while (iter_limit > step)
    {
        error = static_cast<double>(0.0);
        flt_idx = 0;
        for (std::size_t r = 0; r < matrix_data_ptr_->rows(); ++r, flt_idx += cols)
        {
            sigma_1 = sigma_2 = static_cast<double>(0.0);
            diag = matrix_data_ptr_->diagonal_at_row(r);
            cols = matrix_data_ptr_->non_zero_column_size(r);
            for (std::size_t c = flt_idx; c < cols + flt_idx; ++c)
            {
                mat_elm = matrix_data_ptr_->at(c);
                mat_r = std::get<0>(mat_elm);
                mat_c = std::get<1>(mat_elm);
                mat_val = std::get<2>(mat_elm);
                if (mat_c < mat_r)
                    sigma_1 += mat_val * x_new[mat_c];
                if (mat_c > mat_r)
                    sigma_2 += mat_val * solution[mat_c];
            }
            x_new[r] = (one - omega_) * solution[r] + ((omega_ / diag) * (b_[r] - sigma_1 - sigma_2));
            error += (x_new[r] - solution[r]) * (x_new[r] - solution[r]);
        }

        if (error <= tol)
            break;
        solution = x_new;
        step++;
    }
}

void core_sor_solver::solve(container_t &solution)
{
    LSS_ASSERT(solution.size() == system_size_, "Incorrect size of solution container");
    kernel(solution);
}

container_t const core_sor_solver::solve()
{
    container_t solution(system_size_);
    kernel(solution);
    return solution;
}

} // namespace lss_core_sor_solver
