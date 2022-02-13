#include "splitting_config.hpp"

namespace lss
{

splitting_config_builder::splitting_config_builder()
{
}

splitting_config_builder &splitting_config_builder::splitting(splitting_method splitting_method)
{
    splitting_method_ = splitting_method;
    return *this;
}

splitting_config_builder &splitting_config_builder::weighting_value(double value)
{
    weighting_value_ = value;
    return *this;
}

splitting_config_ptr splitting_config_builder::build()
{
    return std::make_shared<splitting_config>(splitting_method_, weighting_value_);
}

} // namespace lss
