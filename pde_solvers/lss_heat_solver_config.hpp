/**

    @file      lss_heat_solver_config.hpp
    @brief     Solver configurations for heat problems
    @details   ~
    @author    Michal Sara
    @date      13.02.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_LSS_HEAT_SOLVER_CONFIG_HPP_)
#define _LSS_HEAT_SOLVER_CONFIG_HPP_

#include "../common/lss_enumerations.hpp"
#include "../common/lss_macros.hpp"
#include "../common/lss_utility.hpp"
#include "lss_implicit_pde_scheme.hpp"
#include "lss_pde_solver_config.hpp"

namespace lss_pde_solvers
{

using lss_enumerations::dimension_enum;
using lss_enumerations::explicit_pde_schemes_enum;
using lss_enumerations::factorization_enum;
using lss_enumerations::memory_space_enum;
using lss_enumerations::traverse_direction_enum;
using lss_enumerations::tridiagonal_method_enum;
using lss_utility::sptr_t;

/**
    heat_implicit_solver_config structure
 */
struct heat_implicit_solver_config : public pde_implicit_solver_config
{
  private:
    double implicit_pde_scheme_value_;

    explicit heat_implicit_solver_config() = delete;

  public:
    /**
        @brief heat_implicit_solver_config object constructor
        @param memory_space
        @param traverse_direction
        @param tridiagonal_method
        @param tridiagonal_factorization
        @param implicit_pde_scheme
    **/
    LSS_API explicit heat_implicit_solver_config(memory_space_enum const &memory_space,
                                                 traverse_direction_enum const &traverse_direction,
                                                 tridiagonal_method_enum const &tridiagonal_method,
                                                 factorization_enum const &tridiagonal_factorization,
                                                 implicit_pde_scheme_ptr const &implicit_pde_scheme);
    LSS_API ~heat_implicit_solver_config();

    /**
        @brief Implicit scheme
        @retval value for implicit scheme
    **/
    LSS_API double implicit_pde_scheme_value() const;
};

/**
    heat_explicit_solver_config structure
 */
struct heat_explicit_solver_config : public pde_explicit_solver_config
{
  private:
    explicit_pde_schemes_enum explicit_pde_scheme_;

    explicit heat_explicit_solver_config() = delete;

  public:
    /**
        @brief heat_explicit_solver_config object constructor
        @param memory_space
        @param traverse_direction
        @param explicit_pde_scheme
    **/
    LSS_API explicit heat_explicit_solver_config(memory_space_enum const &memory_space,
                                                 traverse_direction_enum const &traverse_direction,
                                                 explicit_pde_schemes_enum const &explicit_pde_scheme);
    LSS_API ~heat_explicit_solver_config();

    /**
        @brief  Explicit PDE scheme
        @retval enum for explicit PDE scheme
    **/
    LSS_API explicit_pde_schemes_enum explicit_pde_scheme() const;
};

using heat_implicit_solver_config_ptr = sptr_t<heat_implicit_solver_config>;

using heat_explicit_solver_config_ptr = sptr_t<heat_explicit_solver_config>;

namespace default_heat_solver_configs
{
// =================================================
// ===== some default implicit solver configs ======
// =================================================
// forward stepping:

namespace
{

/**
    @brief  Quickly build implicit config with value for implicit PDE scheme
    @param  memory
    @param  traverse_direction
    @param  tridiagonal_method
    @param  factorization
    @param  value
    @retval
**/
static const heat_implicit_solver_config_ptr build_implicit_config(memory_space_enum memory,
                                                                   traverse_direction_enum traverse_direction,
                                                                   tridiagonal_method_enum tridiagonal_method,
                                                                   factorization_enum factorization, double value)
{
    return std::make_shared<heat_implicit_solver_config>(memory, traverse_direction, tridiagonal_method, factorization,
                                                         std::make_shared<implicit_pde_scheme>(value));
}

/**
    @brief  Quickly build implicit config with implicit enum pde scheme
    @param  memory
    @param  traverse_direction
    @param  tridiagonal_method
    @param  factorization
    @param  scheme
    @retval
**/
static const heat_implicit_solver_config_ptr build_implicit_config(memory_space_enum memory,
                                                                   traverse_direction_enum traverse_direction,
                                                                   tridiagonal_method_enum tridiagonal_method,
                                                                   factorization_enum factorization,
                                                                   implicit_pde_schemes_enum scheme)
{
    if (scheme == implicit_pde_schemes_enum::Euler)
        return std::make_shared<heat_implicit_solver_config>(
            memory, traverse_direction, tridiagonal_method, factorization,
            std::make_shared<implicit_pde_scheme>(implicit_pde_schemes_enum::Euler));
    else
        return std::make_shared<heat_implicit_solver_config>(
            memory, traverse_direction, tridiagonal_method, factorization,
            std::make_shared<implicit_pde_scheme>(implicit_pde_schemes_enum::CrankNicolson));
}

} // namespace

static auto dev_fwd_cusolver_qr_euler_solver_config_ptr = build_implicit_config(
    memory_space_enum::Device, traverse_direction_enum::Forward, tridiagonal_method_enum::CUDASolver,
    factorization_enum::QRMethod, implicit_pde_schemes_enum::Euler);

static auto dev_fwd_cusolver_qr_cn_solver_config_ptr = build_implicit_config(
    memory_space_enum::Device, traverse_direction_enum::Forward, tridiagonal_method_enum::CUDASolver,
    factorization_enum::QRMethod, implicit_pde_schemes_enum::CrankNicolson);

static auto dev_fwd_cusolver_lu_euler_solver_config_ptr = build_implicit_config(
    memory_space_enum::Device, traverse_direction_enum::Forward, tridiagonal_method_enum::CUDASolver,
    factorization_enum::LUMethod, implicit_pde_schemes_enum::Euler);

static auto dev_fwd_cusolver_lu_cn_solver_config_ptr = build_implicit_config(
    memory_space_enum::Device, traverse_direction_enum::Forward, tridiagonal_method_enum::CUDASolver,
    factorization_enum::LUMethod, implicit_pde_schemes_enum::CrankNicolson);

static auto dev_fwd_cusolver_euler_solver_config_ptr = build_implicit_config(
    memory_space_enum::Device, traverse_direction_enum::Forward, tridiagonal_method_enum::CUDASolver,
    factorization_enum::None, implicit_pde_schemes_enum::Euler);

static auto host_fwd_cusolver_qr_euler_solver_config_ptr = build_implicit_config(
    memory_space_enum::Host, traverse_direction_enum::Forward, tridiagonal_method_enum::CUDASolver,
    factorization_enum::QRMethod, implicit_pde_schemes_enum::Euler);

static auto host_fwd_cusolver_qr_cn_solver_config_ptr = build_implicit_config(
    memory_space_enum::Host, traverse_direction_enum::Forward, tridiagonal_method_enum::CUDASolver,
    factorization_enum::QRMethod, implicit_pde_schemes_enum::CrankNicolson);

static auto host_fwd_cusolver_lu_euler_solver_config_ptr = build_implicit_config(
    memory_space_enum::Host, traverse_direction_enum::Forward, tridiagonal_method_enum::CUDASolver,
    factorization_enum::LUMethod, implicit_pde_schemes_enum::Euler);

static auto host_fwd_cusolver_lu_cn_solver_config_ptr = build_implicit_config(
    memory_space_enum::Host, traverse_direction_enum::Forward, tridiagonal_method_enum::CUDASolver,
    factorization_enum::LUMethod, implicit_pde_schemes_enum::CrankNicolson);

static auto host_fwd_cusolver_euler_solver_config_ptr = build_implicit_config(
    memory_space_enum::Host, traverse_direction_enum::Forward, tridiagonal_method_enum::CUDASolver,
    factorization_enum::None, implicit_pde_schemes_enum::Euler);

static auto dev_fwd_sorsolver_euler_solver_config_ptr = build_implicit_config(
    memory_space_enum::Device, traverse_direction_enum::Forward, tridiagonal_method_enum::SORSolver,
    factorization_enum::None, implicit_pde_schemes_enum::Euler);

static auto dev_fwd_sorsolver_cn_solver_config_ptr = build_implicit_config(
    memory_space_enum::Device, traverse_direction_enum::Forward, tridiagonal_method_enum::SORSolver,
    factorization_enum::None, implicit_pde_schemes_enum::CrankNicolson);

static auto host_fwd_sorsolver_euler_solver_config_ptr =
    build_implicit_config(memory_space_enum::Host, traverse_direction_enum::Forward, tridiagonal_method_enum::SORSolver,
                          factorization_enum::None, implicit_pde_schemes_enum::Euler);

static auto host_fwd_sorsolver_cn_solver_config_ptr =
    build_implicit_config(memory_space_enum::Host, traverse_direction_enum::Forward, tridiagonal_method_enum::SORSolver,
                          factorization_enum::None, implicit_pde_schemes_enum::CrankNicolson);

static auto host_fwd_dssolver_euler_solver_config_ptr = build_implicit_config(
    memory_space_enum::Host, traverse_direction_enum::Forward, tridiagonal_method_enum::DoubleSweepSolver,
    factorization_enum::None, implicit_pde_schemes_enum::Euler);

static auto host_fwd_dssolver_cn_solver_config_ptr = build_implicit_config(
    memory_space_enum::Host, traverse_direction_enum::Forward, tridiagonal_method_enum::DoubleSweepSolver,
    factorization_enum::None, implicit_pde_schemes_enum::CrankNicolson);

static auto host_fwd_tlusolver_euler_solver_config_ptr = build_implicit_config(
    memory_space_enum::Host, traverse_direction_enum::Forward, tridiagonal_method_enum::ThomasLUSolver,
    factorization_enum::None, implicit_pde_schemes_enum::Euler);

static auto host_fwd_tlusolver_cn_solver_config_ptr = build_implicit_config(
    memory_space_enum::Host, traverse_direction_enum::Forward, tridiagonal_method_enum::ThomasLUSolver,
    factorization_enum::None, implicit_pde_schemes_enum::CrankNicolson);

// backward stepping:

static auto dev_bwd_cusolver_qr_euler_solver_config_ptr = build_implicit_config(
    memory_space_enum::Device, traverse_direction_enum::Backward, tridiagonal_method_enum::CUDASolver,
    factorization_enum::QRMethod, implicit_pde_schemes_enum::Euler);

static auto dev_bwd_cusolver_qr_cn_solver_config_ptr = build_implicit_config(
    memory_space_enum::Device, traverse_direction_enum::Backward, tridiagonal_method_enum::CUDASolver,
    factorization_enum::QRMethod, implicit_pde_schemes_enum::CrankNicolson);

static auto dev_bwd_cusolver_lu_euler_solver_config_ptr = build_implicit_config(
    memory_space_enum::Device, traverse_direction_enum::Backward, tridiagonal_method_enum::CUDASolver,
    factorization_enum::LUMethod, implicit_pde_schemes_enum::Euler);

static auto dev_bwd_cusolver_lu_cn_solver_config_ptr = build_implicit_config(
    memory_space_enum::Device, traverse_direction_enum::Backward, tridiagonal_method_enum::CUDASolver,
    factorization_enum::LUMethod, implicit_pde_schemes_enum::CrankNicolson);

static auto dev_bwd_cusolver_euler_solver_config_ptr = build_implicit_config(
    memory_space_enum::Device, traverse_direction_enum::Backward, tridiagonal_method_enum::CUDASolver,
    factorization_enum::None, implicit_pde_schemes_enum::Euler);

static auto host_bwd_cusolver_qr_euler_solver_config_ptr = build_implicit_config(
    memory_space_enum::Host, traverse_direction_enum::Backward, tridiagonal_method_enum::CUDASolver,
    factorization_enum::QRMethod, implicit_pde_schemes_enum::Euler);

static auto host_bwd_cusolver_qr_cn_solver_config_ptr = build_implicit_config(
    memory_space_enum::Host, traverse_direction_enum::Backward, tridiagonal_method_enum::CUDASolver,
    factorization_enum::QRMethod, implicit_pde_schemes_enum::CrankNicolson);

static auto host_bwd_cusolver_lu_euler_solver_config_ptr = build_implicit_config(
    memory_space_enum::Host, traverse_direction_enum::Backward, tridiagonal_method_enum::CUDASolver,
    factorization_enum::LUMethod, implicit_pde_schemes_enum::Euler);

static auto host_bwd_cusolver_lu_cn_solver_config_ptr = build_implicit_config(
    memory_space_enum::Host, traverse_direction_enum::Backward, tridiagonal_method_enum::CUDASolver,
    factorization_enum::LUMethod, implicit_pde_schemes_enum::CrankNicolson);

static auto host_bwd_cusolver_euler_solver_config_ptr = build_implicit_config(
    memory_space_enum::Host, traverse_direction_enum::Backward, tridiagonal_method_enum::CUDASolver,
    factorization_enum::None, implicit_pde_schemes_enum::Euler);

static auto dev_bwd_sorsolver_euler_solver_config_ptr = build_implicit_config(
    memory_space_enum::Host, traverse_direction_enum::Backward, tridiagonal_method_enum::SORSolver,
    factorization_enum::None, implicit_pde_schemes_enum::Euler);

static auto dev_bwd_sorsolver_cn_solver_config_ptr = build_implicit_config(
    memory_space_enum::Device, traverse_direction_enum::Backward, tridiagonal_method_enum::SORSolver,
    factorization_enum::None, implicit_pde_schemes_enum::CrankNicolson);

static auto host_bwd_sorsolver_euler_solver_config_ptr = build_implicit_config(
    memory_space_enum::Host, traverse_direction_enum::Backward, tridiagonal_method_enum::SORSolver,
    factorization_enum::None, implicit_pde_schemes_enum::Euler);

static auto host_bwd_sorsolver_cn_solver_config_ptr = build_implicit_config(
    memory_space_enum::Host, traverse_direction_enum::Backward, tridiagonal_method_enum::SORSolver,
    factorization_enum::None, implicit_pde_schemes_enum::CrankNicolson);

static auto host_bwd_dssolver_euler_solver_config_ptr = build_implicit_config(
    memory_space_enum::Host, traverse_direction_enum::Backward, tridiagonal_method_enum::DoubleSweepSolver,
    factorization_enum::None, implicit_pde_schemes_enum::Euler);

static auto host_bwd_dssolver_cn_solver_config_ptr = build_implicit_config(
    memory_space_enum::Host, traverse_direction_enum::Backward, tridiagonal_method_enum::DoubleSweepSolver,
    factorization_enum::None, implicit_pde_schemes_enum::CrankNicolson);

static auto host_bwd_tlusolver_euler_solver_config_ptr = build_implicit_config(
    memory_space_enum::Host, traverse_direction_enum::Backward, tridiagonal_method_enum::ThomasLUSolver,
    factorization_enum::None, implicit_pde_schemes_enum::Euler);

static auto host_bwd_tlusolver_cn_solver_config_ptr = build_implicit_config(
    memory_space_enum::Host, traverse_direction_enum::Backward, tridiagonal_method_enum::ThomasLUSolver,
    factorization_enum::None, implicit_pde_schemes_enum::CrankNicolson);

// =================================================
// ===== some default explicit solver configs ======
// =================================================
static auto dev_expl_fwd_euler_solver_config_ptr = std::make_shared<heat_explicit_solver_config>(
    memory_space_enum::Device, traverse_direction_enum::Forward, explicit_pde_schemes_enum::Euler);

static auto host_expl_fwd_euler_solver_config_ptr = std::make_shared<heat_explicit_solver_config>(
    memory_space_enum::Host, traverse_direction_enum::Forward, explicit_pde_schemes_enum::Euler);

static auto dev_expl_bwd_euler_solver_config_ptr = std::make_shared<heat_explicit_solver_config>(
    memory_space_enum::Device, traverse_direction_enum::Backward, explicit_pde_schemes_enum::Euler);

static auto host_expl_bwd_euler_solver_config_ptr = std::make_shared<heat_explicit_solver_config>(
    memory_space_enum::Host, traverse_direction_enum::Backward, explicit_pde_schemes_enum::Euler);

static auto host_expl_fwd_bc_solver_config_ptr = std::make_shared<heat_explicit_solver_config>(
    memory_space_enum::Host, traverse_direction_enum::Forward, explicit_pde_schemes_enum::ADEBarakatClark);

static auto host_expl_bwd_bc_solver_config_ptr = std::make_shared<heat_explicit_solver_config>(
    memory_space_enum::Host, traverse_direction_enum::Backward, explicit_pde_schemes_enum::ADEBarakatClark);

static auto host_expl_fwd_s_solver_config_ptr = std::make_shared<heat_explicit_solver_config>(
    memory_space_enum::Host, traverse_direction_enum::Forward, explicit_pde_schemes_enum::ADESaulyev);

static auto host_expl_bwd_s_solver_config_ptr = std::make_shared<heat_explicit_solver_config>(
    memory_space_enum::Host, traverse_direction_enum::Backward, explicit_pde_schemes_enum::ADESaulyev);
} // namespace default_heat_solver_configs
} // namespace lss_pde_solvers

#endif ///_LSS_HEAT_SOLVER_CONFIG_HPP_
