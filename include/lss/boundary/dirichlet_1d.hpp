#if !defined(_DIRICHLET_1D_HPP_)
#define _DIRICHLET_1D_HPP_

#include <functional>

#include "../../../boundaries/lss_boundary.hpp"
#include "../../../boundaries/lss_dirichlet_boundary.hpp"
#include "../../../common/lss_macros.hpp"

namespace lss
{

using dirichlet_1d_ptr = lss_boundary::dirichlet_boundary_1d_ptr;
using dirichlet_1d = lss_boundary::dirichlet_boundary_1d;

struct dirichlet_1d_builder
{
  private:
    double const_value_;
    std::function<double(double)> fun_value_{nullptr};

  public:
    LSS_API explicit dirichlet_1d_builder();

    LSS_API dirichlet_1d_builder &value(const std::function<double(double)> &value);

    LSS_API dirichlet_1d_builder &value(double value);

    LSS_API dirichlet_1d_ptr build();
};
} // namespace lss

#endif ///_DIRICHLET_1D_HPP_
