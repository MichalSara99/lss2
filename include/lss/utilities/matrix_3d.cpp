#include "matrix_3d.hpp"

namespace lss
{

matrix_3d_builder::matrix_3d_builder()
{
}

matrix_3d_builder &matrix_3d_builder::rows(std::size_t rows)
{
    rows_ = rows;
    return *this;
}

matrix_3d_builder &matrix_3d_builder::columns(std::size_t columns)
{
    columns_ = columns;
    return *this;
}

matrix_3d_builder &matrix_3d_builder::layers(std::size_t layers)
{
    layers_ = layers;
    return *this;
}

matrix_3d_builder &matrix_3d_builder::default_value(double value)
{
    dvalue_ = value;
    return *this;
}

matrix_3d_ptr matrix_3d_builder::build()
{
    return std::make_shared<matrix_3d>(rows_, columns_, layers_, dvalue_);
}

} // namespace lss
