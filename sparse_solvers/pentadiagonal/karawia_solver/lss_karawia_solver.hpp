/**

    @file      lss_karawia_solver.hpp
    @brief     Pentadiagonal Karawia Solver
    @details   ~
    @author    Michal Sara
    @date      13.02.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_LSS_KARAWIA_SOLVER_HPP_)
#define _LSS_KARAWIA_SOLVER_HPP_

#pragma warning(disable : 4244)

#include <tuple>
#include <type_traits>

#include "../../../boundaries/lss_boundary.hpp"
#include "../../../common/lss_enumerations.hpp"
#include "../../../common/lss_macros.hpp"
#include "../../../common/lss_utility.hpp"
#include "lss_karawia_boundary.hpp"

namespace lss_karawia_solver
{

using lss_boundary::boundary_1d_pair;
using lss_boundary::boundary_2d_pair;
using lss_utility::container_t;
using lss_utility::sptr_t;

class karawia_solver
{

  private:
    std::size_t discretization_size_;
    container_t a_, b_, c_, d_, e_;
    container_t alpha_, beta_, gamma_, mu_, f_;
    karawia_solver_boundary_ptr karawia_boundary_;

    void kernel(boundary_1d_pair const &boundary, boundary_1d_pair const &other_boundary, container_t &solution,
                double time);

    void kernel(boundary_2d_pair const &boundary, boundary_2d_pair const &other_boundary, container_t &solution,
                double time, double space_args);

    void initialize();

    explicit karawia_solver() = delete;
    // bool is_diagonally_dominant() const;

  public:
    explicit karawia_solver(std::size_t discretization_size);

    ~karawia_solver();

    void set_diagonals(container_t lowest_diagonal, container_t lower_diagonal, container_t diagonal,
                       container_t upper_diagonal, container_t uppest_diagonal);

    void set_rhs(container_t const &rhs);

    void solve(boundary_1d_pair const &boundary, boundary_1d_pair const &other_boundary, container_t &solution);

    void solve(boundary_1d_pair const &boundary, boundary_1d_pair const &other_boundary, container_t &solution,
               double at_time);

    void solve(boundary_2d_pair const &boundary, boundary_2d_pair const &other_boundary, container_t &solution,
               double at_time, double space_arg);
};

using karawia_solver_ptr = sptr_t<karawia_solver>;

} // namespace lss_karawia_solver

#endif ///_LSS_KARAWIA_SOLVER_HPP_
