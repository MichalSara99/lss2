/**

    @file      lss_cuda_solver.hpp
    @brief     Tridiagonal CUDA solver
    @details   ~
    @author    Michal Sara
    @date      12.02.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_LSS_CUDA_SOLVER_HPP_)
#define _LSS_CUDA_SOLVER_HPP_

#include "../../../boundaries/lss_boundary.hpp"
#include "../../../common/lss_enumerations.hpp"
#include "../../../common/lss_utility.hpp"
#include "../../core/cuda_solver/lss_core_cuda_solver.hpp"
#include "../lss_tridiagonal_solver.hpp"
#include "lss_cuda_boundary.hpp"

namespace lss_cuda_solver
{

using lss_boundary::boundary_1d_pair;
using lss_boundary::boundary_2d_pair;
using lss_boundary::boundary_3d_pair;
using lss_core_cuda_solver::real_sparse_solver_cuda;
using lss_enumerations::factorization_enum;
using lss_enumerations::memory_space_enum;
using lss_utility::container_t;
using lss_utility::sptr_t;

template <memory_space_enum memory_space> class cuda_solver : public lss_tridiagonal_solver::tridiagonal_solver
{
  private:
    cuda_boundary_ptr cuda_boundary_;

    void kernel(boundary_1d_pair const &boundary, container_t &solution, factorization_enum factorization,
                double time) override;

    void kernel(boundary_2d_pair const &boundary, container_t &solution, factorization_enum factorization, double time,
                double space_arg) override;

    void kernel(boundary_3d_pair const &boundary, container_t &solution, factorization_enum factorization, double time,
                double space_1_arg, double space_2_arg) override;

    void initialize();

    explicit cuda_solver() = delete;

  public:
    explicit cuda_solver(std::size_t discretization_size);

    ~cuda_solver();

    cuda_solver(cuda_solver const &) = delete;
    cuda_solver(cuda_solver &&) = delete;
    cuda_solver &operator=(cuda_solver const &) = delete;
    cuda_solver &operator=(cuda_solver &&) = delete;
};

using device_cuda_solver_ptr = sptr_t<cuda_solver<memory_space_enum::Device>>;
using host_cuda_solver_ptr = sptr_t<cuda_solver<memory_space_enum::Host>>;

} // namespace lss_cuda_solver

#endif ///_LSS_CUDA_SOLVER_HPP_
