#if !defined(_OUTPUT_HPP_)
#define _OUTPUT_HPP_

#include <iostream>
#include <vector>

#include "../../../common/lss_macros.hpp"
#include "../../../common/lss_print.hpp"
#include "../../../common/lss_utility.hpp"
#include "../../../common/lss_xml.hpp"
#include "../../../discretization/lss_discretization_config.hpp"
#include "../../lss/configs/discretization_config_1d.hpp"
#include "../../lss/configs/discretization_config_2d.hpp"
#include "../../lss/configs/grid_config_1d.hpp"
#include "../../lss/configs/grid_config_2d.hpp"
#include "../../lss/utilities/container_2d.hpp"

namespace lss
{

enum class to
{
    xml,
    stream
};

using discretization_config_ptr = lss_discretization::discretization_config_1d_ptr;

LSS_API extern void output_curve(discretization_config_ptr const &discretization_config,
                                 grid_config_1d_ptr const &grid_config_hints, std::vector<double> const &container,
                                 to to, std::ostream &out = std::cout);

LSS_API extern void output_surface(discretization_config_1d_ptr const &discretization_config,
                                   grid_config_1d_ptr const &grid_config_hints, rmatrix_2d const &container, to to,
                                   std::ostream &out = std::cout);

LSS_API extern void output_surface(discretization_config_2d_ptr const &discretization_config,
                                   grid_config_2d_ptr const &grid_config_hints, rmatrix_2d const &container, to to,
                                   std::ostream &out = std::cout);

} // namespace lss

#endif ///_OUTPUT_HPP_
