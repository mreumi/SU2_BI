#pragma once

#include <array>
#include <vector>
#include <cmath>

#include "COutput.hpp"

class CMeshPoint {
public:
    std::vector<su2double> coord{};   // all coords
    std::vector<su2double> values;        // associated values

    // Default constructor e.g. to call: Point pt1;
    CMeshPoint() = default;
    
    // Construct from an existing std::array of coordinates
    CMeshPoint(const std::vector<su2double>& c, const std::vector<su2double>& v = {})
        : coord(c), values(v) {}
    
    // ---- Copy/move operations ----
    CMeshPoint(const CMeshPoint&)                   = default;
    CMeshPoint(CMeshPoint&&) noexcept               = default;
    CMeshPoint& operator=(const CMeshPoint&)        = default;
    CMeshPoint& operator=(CMeshPoint&&) noexcept    = default;

    // Destructor: also default (or omit entirely)
    ~CMeshPoint() = default;

    // ---- Coordinate access helpers ----

    // Access coordinate by axis (0..Dim-1): p[axis]
    su2double& operator[](std::size_t axis) {
        return coord[axis];
    }
    const su2double& operator[](std::size_t axis) const {
        return coord[axis];
    }

    // ---- Distance (nD) ----
    su2double distance(const CMeshPoint& other) const {
        su2double sum = 0;
        for (std::size_t i = 0; i < coord.size(); ++i) {
            su2double d = coord[i] - other.coord[i];
            sum += d * d;
        }
        return std::sqrt(sum);
    }
};
