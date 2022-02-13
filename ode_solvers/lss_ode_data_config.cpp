#include "lss_ode_data_config.hpp"

#include "../common/lss_enumerations.hpp"
#include "../common/lss_macros.hpp"
#include "../common/lss_utility.hpp"

namespace lss_ode_solvers
{
ode_coefficient_data_config::ode_coefficient_data_config(std::function<double(double)> const &a_coefficient,
                                                         std::function<double(double)> const &b_coefficient)
    : a_coeff_{a_coefficient}, b_coeff_{b_coefficient}
{
    initialize();
}

std::function<double(double)> const &ode_coefficient_data_config::a_coefficient() const
{
    return a_coeff_;
}

std::function<double(double)> const &ode_coefficient_data_config::b_coefficient() const
{
    return b_coeff_;
}

void ode_coefficient_data_config::initialize()
{
    LSS_VERIFY(a_coeff_, "a_coefficient must not be null");
    LSS_VERIFY(b_coeff_, "b_coefficient must not be null");
}

ode_nonhom_data_config::ode_nonhom_data_config(std::function<double(double)> const &nonhom_fun)
    : nonhom_fun_{nonhom_fun}
{
    LSS_VERIFY(nonhom_fun_, "nonhom_fun must not be null");
}

std::function<double(double)> const &ode_nonhom_data_config::nonhom_function() const
{
    return nonhom_fun_;
}

void ode_data_config::initialize()
{
    LSS_VERIFY(coefficient_data_cfg_, "coefficient_data_config must not be null");
}

ode_data_config::ode_data_config(ode_coefficient_data_config_ptr const &coefficient_data_config,
                                 ode_nonhom_data_config_ptr const &nonhom_data_config)
    : coefficient_data_cfg_{coefficient_data_config}, nonhom_data_cfg_{nonhom_data_config}
{
    initialize();
}

ode_data_config ::~ode_data_config()
{
}

ode_nonhom_data_config_ptr const &ode_data_config::nonhom_data_config() const
{
    return nonhom_data_cfg_;
}

std::function<double(double)> const &ode_data_config::a_coefficient() const
{
    return coefficient_data_cfg_->a_coefficient();
}

std::function<double(double)> const &ode_data_config::b_coefficient() const
{
    return coefficient_data_cfg_->b_coefficient();
}
} // namespace lss_ode_solvers
