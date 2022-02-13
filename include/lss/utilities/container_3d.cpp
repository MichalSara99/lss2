#include "container_3d.hpp"

namespace lss
{

rmatrix_3d_builder::rmatrix_3d_builder()
{
}

rmatrix_3d_builder &rmatrix_3d_builder::rows(std::size_t rows)
{
    rows_ = rows;
    return *this;
}

rmatrix_3d_builder &rmatrix_3d_builder::columns(std::size_t columns)
{
    columns_ = columns;
    return *this;
}

rmatrix_3d_builder &rmatrix_3d_builder::layers(std::size_t layers)
{
    layers_ = layers;
    return *this;
}

rmatrix_3d_builder &rmatrix_3d_builder::default_value(double value)
{
    dvalue_ = value;
    return *this;
}

rmatrix_3d_ptr rmatrix_3d_builder::build()
{
    return std::make_shared<rmatrix_3d>(rows_, columns_, layers_, dvalue_);
}

cmatrix_3d_builder::cmatrix_3d_builder()
{
}

cmatrix_3d_builder &cmatrix_3d_builder::rows(std::size_t rows)
{
    rows_ = rows;
    return *this;
}

cmatrix_3d_builder &cmatrix_3d_builder::columns(std::size_t columns)
{
    columns_ = columns;
    return *this;
}

cmatrix_3d_builder &cmatrix_3d_builder::layers(std::size_t layers)
{
    layers_ = layers;
    return *this;
}

cmatrix_3d_builder &cmatrix_3d_builder::default_value(double value)
{
    dvalue_ = value;
    return *this;
}

cmatrix_3d_ptr cmatrix_3d_builder::build()
{
    return std::make_shared<cmatrix_3d>(rows_, columns_, layers_, dvalue_);
}

} // namespace lss
