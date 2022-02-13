#include "lss_wave_data_config.hpp"
#include "../common/lss_macros.hpp"

namespace lss_pde_solvers
{

void wave_coefficient_data_config_1d::initialize()
{
    LSS_VERIFY(a_coeff_, "a_coefficient must not be null");
    LSS_VERIFY(b_coeff_, "b_coefficient must not be null");
    LSS_VERIFY(c_coeff_, "c_coefficient must not be null");
    LSS_VERIFY(d_coeff_, "c_coefficient must not be null");
}

wave_coefficient_data_config_1d ::wave_coefficient_data_config_1d(
    std::function<double(double, double)> const &a_coefficient,
    std::function<double(double, double)> const &b_coefficient,
    std::function<double(double, double)> const &c_coefficient,
    std::function<double(double, double)> const &d_coefficient)
    : a_coeff_{a_coefficient}, b_coeff_{b_coefficient}, c_coeff_{c_coefficient}, d_coeff_{d_coefficient}
{
    initialize();
}

std::function<double(double, double)> const &wave_coefficient_data_config_1d::a_coefficient() const
{
    return a_coeff_;
}

std::function<double(double, double)> const &wave_coefficient_data_config_1d::b_coefficient() const
{
    return b_coeff_;
}

std::function<double(double, double)> const &wave_coefficient_data_config_1d::c_coefficient() const
{
    return c_coeff_;
}

std::function<double(double, double)> const &wave_coefficient_data_config_1d::d_coefficient() const
{
    return d_coeff_;
}

wave_initial_data_config_1d::wave_initial_data_config_1d(std::function<double(double)> const &first_initial_condition,
                                                         std::function<double(double)> const &second_initial_condition)
    : first_initial_condition_{first_initial_condition}, second_initial_condition_{second_initial_condition}
{
    LSS_VERIFY(first_initial_condition_, "first_initial_condition must not be null");
    LSS_VERIFY(second_initial_condition_, "second_initial_condition must not be null");
}

std::function<double(double)> const &wave_initial_data_config_1d::first_initial_condition() const
{
    return first_initial_condition_;
}

std::function<double(double)> const &wave_initial_data_config_1d::second_initial_condition() const
{
    return second_initial_condition_;
}

wave_source_data_config_1d::wave_source_data_config_1d(std::function<double(double, double)> const &wave_source)
    : wave_source_{wave_source}
{
    LSS_VERIFY(wave_source_, "wave_source must not be null");
}

std::function<double(double, double)> const &wave_source_data_config_1d::wave_source() const
{
    return wave_source_;
}

void wave_data_config_1d::initialize()
{
    LSS_VERIFY(coefficient_data_cfg_, "coefficient_data_config must not be null");
    LSS_VERIFY(initial_data_cfg_, "initial_data_config must not be null");
}

wave_data_config_1d::wave_data_config_1d(wave_coefficient_data_config_1d_ptr const &coefficient_data_config,
                                         wave_initial_data_config_1d_ptr const &initial_data_config,
                                         wave_source_data_config_1d_ptr const &source_data_config)
    : coefficient_data_cfg_{coefficient_data_config}, initial_data_cfg_{initial_data_config}, source_data_cfg_{
                                                                                                  source_data_config}
{
    initialize();
}

wave_data_config_1d::~wave_data_config_1d()
{
}

wave_source_data_config_1d_ptr const &wave_data_config_1d::source_data_config() const
{
    return source_data_cfg_;
}

std::function<double(double)> const &wave_data_config_1d::first_initial_condition() const
{
    return initial_data_cfg_->first_initial_condition();
}

std::function<double(double)> const &wave_data_config_1d::second_initial_condition() const
{
    return initial_data_cfg_->second_initial_condition();
}

std::function<double(double, double)> const &wave_data_config_1d::a_coefficient() const
{
    return coefficient_data_cfg_->a_coefficient();
}

std::function<double(double, double)> const &wave_data_config_1d::b_coefficient() const
{
    return coefficient_data_cfg_->b_coefficient();
}

std::function<double(double, double)> const &wave_data_config_1d::c_coefficient() const
{
    return coefficient_data_cfg_->c_coefficient();
}

std::function<double(double, double)> const &wave_data_config_1d::d_coefficient() const
{
    return coefficient_data_cfg_->d_coefficient();
}
} // namespace lss_pde_solvers
