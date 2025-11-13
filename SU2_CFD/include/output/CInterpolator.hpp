#pragma once
#include <string>
#include <vector>
#include "CDatapointcloud.hpp"   // Point, Vec, Neighbor, Pointcloud

// Stateless helper stays available if you like calling it directly.
CDatapoint interpolate_IDW(const CDatapointcloud& pc,
                      const CDatapoint& query_point,
                      int k,
                      su2double power);

class Interpolator {
public:
    // --- Data ---
    int k{3};                   // number of neighbors to consider
    su2double power{2.0};          // power parameter for IDW
    std::string method{"IDW"};  // interpolation method

    // --- Constructors ---
    Interpolator() = default;
    Interpolator(int k_, su2double power_, std::string method_ = "IDW") : k(k_), power(power_), method(method_) {}

    // --- Methods ---
    // Single query
    CDatapoint interpolate(const CDatapointcloud& pc,
                      const CDatapoint& query_point,
                      const std::string& method_in = "") const;
    
    // Many queries
    CDatapointcloud interpolate(const CDatapointcloud& pc,
                                   const CDatapointcloud& queries,
                                   const std::string& method_in = "") const;

};