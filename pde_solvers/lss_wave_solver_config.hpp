/**

    @file      lss_wave_solver_config.hpp
    @brief     Solver configuration objects for wave problems
    @details   ~
    @author    Michal Sara
    @date      13.02.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_LSS_WAVE_SOLVER_CONFIG_HPP_)
#define _LSS_WAVE_SOLVER_CONFIG_HPP_

#include <map>

#include "../common/lss_enumerations.hpp"
#include "../common/lss_utility.hpp"
#include "lss_pde_solver_config.hpp"

namespace lss_pde_solvers
{

using lss_enumerations::explicit_pde_schemes_enum;
using lss_enumerations::factorization_enum;
using lss_enumerations::implicit_pde_schemes_enum;
using lss_enumerations::memory_space_enum;
using lss_enumerations::traverse_direction_enum;
using lss_enumerations::tridiagonal_method_enum;
using lss_utility::sptr_t;

/**
    wave_implicit_solver_config structure
 */
struct wave_implicit_solver_config : public pde_implicit_solver_config
{
    explicit wave_implicit_solver_config() = delete;

  public:
    explicit wave_implicit_solver_config(memory_space_enum const &memory_space,
                                         traverse_direction_enum const &traverse_direction,
                                         tridiagonal_method_enum const &tridiagonal_method,
                                         factorization_enum const &tridiagonal_factorization);
    ~wave_implicit_solver_config();
};

/**
    wave_explicit_solver_config structure
 */
struct wave_explicit_solver_config : public pde_explicit_solver_config
{

    explicit wave_explicit_solver_config() = delete;

  public:
    explicit wave_explicit_solver_config(memory_space_enum const &memory_space,
                                         traverse_direction_enum const &traverse_direction);
    ~wave_explicit_solver_config();
};

using wave_implicit_solver_config_ptr = sptr_t<wave_implicit_solver_config>;

using wave_explicit_solver_config_ptr = sptr_t<wave_explicit_solver_config>;

namespace default_wave_solver_configs
{
// =================================================
// ===== some default implicit solver configs ======
// =================================================
// forward stepping:

static auto dev_fwd_cusolver_qr_solver_config_ptr =
    std::make_shared<wave_implicit_solver_config>(memory_space_enum::Device, traverse_direction_enum::Forward,
                                                  tridiagonal_method_enum::CUDASolver, factorization_enum::QRMethod);

static auto dev_fwd_sorsolver_solver_config_ptr =
    std::make_shared<wave_implicit_solver_config>(memory_space_enum::Device, traverse_direction_enum::Forward,
                                                  tridiagonal_method_enum::SORSolver, factorization_enum::None);

static auto host_fwd_sorsolver_solver_config_ptr =
    std::make_shared<wave_implicit_solver_config>(memory_space_enum::Host, traverse_direction_enum::Forward,
                                                  tridiagonal_method_enum::SORSolver, factorization_enum::None);

static auto host_fwd_dssolver_solver_config_ptr =
    std::make_shared<wave_implicit_solver_config>(memory_space_enum::Host, traverse_direction_enum::Forward,
                                                  tridiagonal_method_enum::DoubleSweepSolver, factorization_enum::None);

static auto host_fwd_tlusolver_solver_config_ptr =
    std::make_shared<wave_implicit_solver_config>(memory_space_enum::Host, traverse_direction_enum::Forward,
                                                  tridiagonal_method_enum::ThomasLUSolver, factorization_enum::None);

// backward stepping:

static auto dev_bwd_cusolver_qr_solver_config_ptr =
    std::make_shared<wave_implicit_solver_config>(memory_space_enum::Device, traverse_direction_enum::Backward,
                                                  tridiagonal_method_enum::CUDASolver, factorization_enum::QRMethod);

static auto dev_bwd_cusolver_lu_solver_config_ptr =
    std::make_shared<wave_implicit_solver_config>(memory_space_enum::Device, traverse_direction_enum::Backward,
                                                  tridiagonal_method_enum::CUDASolver, factorization_enum::LUMethod);

static auto dev_bwd_sorsolver_solver_config_ptr =
    std::make_shared<wave_implicit_solver_config>(memory_space_enum::Device, traverse_direction_enum::Backward,
                                                  tridiagonal_method_enum::SORSolver, factorization_enum::None);

static auto host_bwd_sorsolver_solver_config_ptr =
    std::make_shared<wave_implicit_solver_config>(memory_space_enum::Host, traverse_direction_enum::Backward,
                                                  tridiagonal_method_enum::SORSolver, factorization_enum::None);

static auto host_bwd_dssolver_solver_config_ptr =
    std::make_shared<wave_implicit_solver_config>(memory_space_enum::Host, traverse_direction_enum::Backward,
                                                  tridiagonal_method_enum::DoubleSweepSolver, factorization_enum::None);

static auto host_bwd_tlusolver_solver_config_ptr =
    std::make_shared<wave_implicit_solver_config>(memory_space_enum::Host, traverse_direction_enum::Backward,
                                                  tridiagonal_method_enum::ThomasLUSolver, factorization_enum::None);

// =================================================
// ===== some default explicit solver configs ======
// =================================================
// forward stepping:

static auto dev_expl_fwd_solver_config_ptr =
    std::make_shared<wave_explicit_solver_config>(memory_space_enum::Device, traverse_direction_enum::Forward);

static auto host_expl_fwd_solver_config_ptr =
    std::make_shared<wave_explicit_solver_config>(memory_space_enum::Host, traverse_direction_enum::Forward);

// backward stepping:

static auto dev_expl_bwd_solver_config_ptr =
    std::make_shared<wave_explicit_solver_config>(memory_space_enum::Device, traverse_direction_enum::Backward);

static auto host_expl_bwd_solver_config_ptr =
    std::make_shared<wave_explicit_solver_config>(memory_space_enum::Host, traverse_direction_enum::Backward);

} // namespace default_wave_solver_configs

} // namespace lss_pde_solvers

#endif ///_LSS_WAVE_SOLVER_CONFIG_HPP_
