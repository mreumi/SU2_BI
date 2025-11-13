#pragma once

#include <vector>
#include <cmath>
#include <stdexcept>

#include "COutput.hpp"

// Type aliases for convenience
using Vec  = std::vector<su2double>;
using Mat  = std::vector<Vec>;

class CDatapoint {
public:
    // --- Data ---
    Vec coords;
    Vec values;

    // --- Constructors ---
    CDatapoint();
    explicit CDatapoint(const Vec& c);
    CDatapoint(const Vec& c, const Vec& v);

    // --- Methods ---
    su2double distance_to(const CDatapoint& other) const;
};
