#include "../../include/output/CDatapoint.hpp"    // include the header

#include <cmath>      // for std::sqrt
#include <stdexcept>  // for std::runtime_error

// --- Constructors ---
CDatapoint::CDatapoint() = default;

CDatapoint::CDatapoint(const Vec& c)
    : coords(c)
{}

CDatapoint::CDatapoint(const Vec& c, const Vec& v)
    : coords(c), values(v)
{}

// --- Methods ---
su2double CDatapoint::distance_to(const CDatapoint& other) const {
    if (coords.size() != other.coords.size()) {
        throw std::runtime_error("Point::distance: dimension mismatch");
    }

    su2double sum = 0.0;
    for (std::size_t i = 0; i < coords.size(); ++i) {
        su2double d = coords[i] - other.coords[i];
        sum += d * d;
    }
    return std::sqrt(sum);
}
