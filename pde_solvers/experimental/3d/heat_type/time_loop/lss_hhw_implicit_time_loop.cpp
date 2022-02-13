#include "lss_hhw_implicit_time_loop.hpp"

#include "../../../../../boundaries/lss_dirichlet_boundary.hpp"

namespace lss_pde_solvers
{

using lss_boundary::dirichlet_boundary_3d;

namespace three_dimensional
{

boundary_3d_pair hhw_implicit_boundaries::get_y(grid_config_2d_ptr const &grid_config_xz,
                                                container_3d<by_enum::RowPlane> const &next_solution)
{
    auto const lci = next_solution.columns() - 1;
    auto lower = [=](double t, double x, double z) -> double {
        const std::size_t i = grid_config_xz->index_of_1(x);
        const std::size_t k = grid_config_xz->index_of_2(z);
        return next_solution(i, 0, k);
    };
    auto upper = [=](double t, double x, double z) -> double {
        const std::size_t i = grid_config_xz->index_of_1(x);
        const std::size_t k = grid_config_xz->index_of_2(z);
        return next_solution(i, lci, k);
    };
    auto y_low = std::make_shared<dirichlet_boundary_3d>(lower);
    auto y_high = std::make_shared<dirichlet_boundary_3d>(upper);
    return std::make_pair(y_low, y_high);
}

boundary_3d_pair hhw_implicit_boundaries::get_intermed_x(grid_config_2d_ptr const &grid_config_yz,
                                                         container_3d<by_enum::RowPlane> const &next_solution)
{
    const std::size_t lri = next_solution.rows() - 1;
    auto lower = [=](double t, double y, double z) -> double {
        const std::size_t j = grid_config_yz->index_of_1(y);
        const std::size_t k = grid_config_yz->index_of_2(z);
        return next_solution(0, j, k);
    };
    auto upper = [=](double t, double y, double z) -> double {
        const std::size_t j = grid_config_yz->index_of_1(y);
        const std::size_t k = grid_config_yz->index_of_2(z);
        return next_solution(lri, j, k);
    };
    auto x_low = std::make_shared<dirichlet_boundary_3d>(lower);
    auto x_high = std::make_shared<dirichlet_boundary_3d>(upper);
    return std::make_pair(x_low, x_high);
}

boundary_3d_pair hhw_implicit_boundaries::get_intermed_z(grid_config_2d_ptr const &grid_config_xy,
                                                         container_3d<by_enum::RowPlane> const &next_solution)
{
    const std::size_t lli = next_solution.layers() - 1;
    auto lower = [=](double t, double x, double y) -> double {
        const std::size_t i = grid_config_xy->index_of_1(x);
        const std::size_t j = grid_config_xy->index_of_2(y);
        return next_solution(i, j, 0);
    };
    auto upper = [=](double t, double x, double y) -> double {
        const std::size_t i = grid_config_xy->index_of_1(x);
        const std::size_t j = grid_config_xy->index_of_2(y);
        return next_solution(i, j, lli);
    };
    auto z_low = std::make_shared<dirichlet_boundary_3d>(lower);
    auto z_high = std::make_shared<dirichlet_boundary_3d>(upper);
    return std::make_pair(z_low, z_high);
}

void hhw_implicit_time_loop::run(heat_splitting_method_3d_ptr const &solver_ptr,
                                 hhw_explicit_boundary_solver_ptr const &boundary_solver_ptr,
                                 boundary_3d_pair const &x_boundary_pair, boundary_3d_ptr const &y_upper_boundary_ptr,
                                 boundary_3d_pair const &z_boundary_pair, grid_config_3d_ptr const &grid_config,
                                 range_ptr const &time_range, std::size_t const &last_time_idx, double const time_step,
                                 traverse_direction_enum const &traverse_dir,
                                 container_3d<by_enum::RowPlane> &prev_solution,
                                 container_3d<by_enum::RowPlane> &next_solution)
{

    const double start_time = time_range->lower();
    const double end_time = time_range->upper();
    const double k = time_step;
    boundary_3d_pair y_boundary_pair;
    boundary_3d_pair x_inter_boundary_pair;
    boundary_3d_pair z_inter_boundary_pair;

    double time{};
    std::size_t time_idx{};

    if (traverse_dir == traverse_direction_enum::Forward)
    {
        time = start_time + k;
        time_idx = 1;
        while (time_idx <= last_time_idx)
        {
            boundary_solver_ptr->solve(prev_solution, x_boundary_pair, y_upper_boundary_ptr, z_boundary_pair, time,
                                       next_solution);
            y_boundary_pair = hhw_implicit_boundaries::get_y(grid_config->grid_13(), next_solution);
            x_inter_boundary_pair = hhw_implicit_boundaries::get_intermed_x(grid_config->grid_23(), prev_solution);
            z_inter_boundary_pair = hhw_implicit_boundaries::get_intermed_z(grid_config->grid_12(), prev_solution);
            solver_ptr->solve(prev_solution, x_inter_boundary_pair, y_boundary_pair, z_inter_boundary_pair, time,
                              next_solution);
            boundary_solver_ptr->solve(prev_solution, x_boundary_pair, z_boundary_pair, time, next_solution);

            prev_solution = next_solution;
            time += k;
            time_idx++;
        }
    }
    else
    {
        time = end_time - k;
        time_idx = last_time_idx;
        do
        {
            time_idx--;
            boundary_solver_ptr->solve(prev_solution, x_boundary_pair, y_upper_boundary_ptr, z_boundary_pair, time,
                                       next_solution);
            y_boundary_pair = hhw_implicit_boundaries::get_y(grid_config->grid_13(), next_solution);
            x_inter_boundary_pair = hhw_implicit_boundaries::get_intermed_x(grid_config->grid_23(), prev_solution);
            z_inter_boundary_pair = hhw_implicit_boundaries::get_intermed_z(grid_config->grid_12(), prev_solution);
            solver_ptr->solve(prev_solution, x_inter_boundary_pair, y_boundary_pair, z_inter_boundary_pair, time,
                              next_solution);
            // std::cout << "Before:\n";
            // for (std::size_t l = 0; l < next_solution.layers(); ++l)
            //{
            //     for (std::size_t r = 0; r < next_solution.rows(); ++r)
            //     {
            //         for (std::size_t c = 0; c < next_solution.columns(); ++c)
            //         {
            //             std::cout << next_solution(r, c, l) << ",";
            //         }
            //         std::cout << "\n";
            //     }
            // }

            boundary_solver_ptr->solve(prev_solution, x_boundary_pair, z_boundary_pair, time, next_solution);
            // std::cout << "After:\n";
            // for (std::size_t l = 0; l < next_solution.layers(); ++l)
            //{
            //     for (std::size_t r = 0; r < next_solution.rows(); ++r)
            //     {
            //         for (std::size_t c = 0; c < next_solution.columns(); ++c)
            //         {
            //             std::cout << next_solution(r, c, l) << ",";
            //         }
            //         std::cout << "\n";
            //     }
            // }

            prev_solution = next_solution;
            time -= k;
        } while (time_idx > 0);
    }
}

} // namespace three_dimensional
} // namespace lss_pde_solvers
