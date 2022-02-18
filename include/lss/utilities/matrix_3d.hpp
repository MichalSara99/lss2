#if !defined(_MATRIX_3D_HPP_)
#define _MATRIX_3D_HPP_

#include <functional>

#include "../../../common/lss_macros.hpp"
#include "../../../common/lss_utility.hpp"
#include "../../../containers/lss_matrix_3d.hpp"

namespace lss
{

using matrix_3d = lss_containers::matrix_3d;
using matrix_3d_ptr = lss_containers::matrix_3d_ptr;

struct matrix_3d_builder
{
  private:
    std::size_t rows_;
    std::size_t columns_;
    std::size_t layers_;
    double dvalue_;

  public:
    LSS_API explicit matrix_3d_builder();

    LSS_API matrix_3d_builder &rows(std::size_t rows);

    LSS_API matrix_3d_builder &columns(std::size_t columns);

    LSS_API matrix_3d_builder &layers(std::size_t layers);

    LSS_API matrix_3d_builder &default_value(double value);

    LSS_API matrix_3d_ptr build();
};

} // namespace lss

#endif ///_MATRIX_3D_HPP_
