#if !defined(_ROBIN_1D_HPP_)
#define _ROBIN_1D_HPP_

#include <functional>

#include "../../../boundaries/lss_boundary.hpp"
#include "../../../boundaries/lss_robin_boundary.hpp"
#include "../../../common/lss_macros.hpp"

namespace lss
{

using robin_1d_ptr = lss_boundary::robin_boundary_1d_ptr;
using robin_1d = lss_boundary::robin_boundary_1d;

struct robin_1d_builder
{
  private:
    double const_linear_value_, const_value_;
    std::function<double(double)> fun_linear_value_;
    std::function<double(double)> fun_value_;

  public:
    LSS_API explicit robin_1d_builder();

    LSS_API robin_1d_builder &values(double linear_value, double value);

    LSS_API robin_1d_builder &values(const std::function<double(double)> &linear_value,
                                     const std::function<double(double)> &value);

    LSS_API robin_1d_ptr build();
};

} // namespace lss

#endif ///_ROBIN_1D_HPP_
