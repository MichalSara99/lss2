/**

    @file      lss_thomas_lu_solver.hpp
    @brief     Tridiagonal Thomas LU Solver
    @details   ~
    @author    Michal Sara
    @date      13.02.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_LSS_THOMAS_LU_SOLVER_HPP_)
#define _LSS_THOMAS_LU_SOLVER_HPP_

#pragma warning(disable : 4244)

#include <tuple>
#include <type_traits>
#include <vector>

#include "../../../boundaries/lss_boundary.hpp"
#include "../../../common/lss_utility.hpp"
#include "../../tridiagonal/lss_tridiagonal_solver.hpp"
#include "lss_thomas_lu_boundary.hpp"

namespace lss_thomas_lu_solver
{

using lss_boundary::boundary_1d_pair;
using lss_boundary::boundary_2d_pair;
using lss_boundary::boundary_3d_pair;
using lss_enumerations::factorization_enum;
using lss_utility::container_t;
using lss_utility::sptr_t;

class thomas_lu_solver : public lss_tridiagonal_solver::tridiagonal_solver
{

  private:
    container_t beta_, gamma_;
    thomas_lu_solver_boundary_ptr tlu_boundary_;

    void kernel(boundary_1d_pair const &boundary, container_t &solution, factorization_enum factorization,
                double time) override;

    void kernel(boundary_2d_pair const &boundary, container_t &solution, factorization_enum factorization, double time,
                double space_arg) override;

    void kernel(boundary_3d_pair const &boundary, container_t &solution, factorization_enum factorization, double time,
                double space_1_arg, double space_2_arg) override;

    void initialize();

    explicit thomas_lu_solver() = delete;
    bool is_diagonally_dominant() const;

  public:
    explicit thomas_lu_solver(std::size_t discretization_size);

    ~thomas_lu_solver();
};

using thomas_lu_solver_ptr = sptr_t<thomas_lu_solver>;

} // namespace lss_thomas_lu_solver

#endif ///_LSS_THOMAS_LU_SOLVER_HPP_
