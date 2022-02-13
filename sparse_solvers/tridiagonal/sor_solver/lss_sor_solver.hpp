/**

    @file      lss_sor_solver.hpp
    @brief     Tridiagonal SOR Solver
    @details   ~
    @author    Michal Sara
    @date      13.02.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_LSS_SOR_SOLVER_HPP_)
#define _LSS_SOR_SOLVER_HPP_

#include "../../../boundaries/lss_boundary.hpp"
#include "../../../common/lss_utility.hpp"
#include "../lss_tridiagonal_solver.hpp"
#include "lss_sor_boundary.hpp"

namespace lss_sor_solver
{

using lss_boundary::boundary_1d_pair;
using lss_boundary::boundary_2d_pair;
using lss_boundary::boundary_3d_pair;
using lss_enumerations::factorization_enum;
using lss_utility::container_t;
using lss_utility::sptr_t;

class sor_solver : public lss_tridiagonal_solver::tridiagonal_solver
{
  private:
    sor_boundary_ptr sor_boundary_;

    void kernel(boundary_1d_pair const &boundary, container_t &solution, factorization_enum factorization,
                double time) override;

    void kernel(boundary_2d_pair const &boundary, container_t &solution, factorization_enum factorization, double time,
                double space_arg) override;

    void kernel(boundary_3d_pair const &boundary, container_t &solution, factorization_enum factorization, double time,
                double space_1_arg, double space_2_arg) override;

    void initialize();

    explicit sor_solver() = delete;

  public:
    explicit sor_solver(std::size_t discretization_size);

    ~sor_solver();

    sor_solver(sor_solver const &) = delete;
    sor_solver(sor_solver &&) = delete;
    sor_solver &operator=(sor_solver const &) = delete;
    sor_solver &operator=(sor_solver &&) = delete;
};

using sor_solver_ptr = sptr_t<sor_solver>;

} // namespace lss_sor_solver

#endif ///_LSS_SOR_SOLVER_HPP_
