// CDatapointcloud.h
#pragma once

#include <vector>
#include <cstddef>   // std::size_t
#include "CDatapoint.hpp"

// Type aliases for convenience
using Vec  = std::vector<su2double>;
using Mat  = std::vector<Vec>;

struct Neighbor {
    CDatapoint point;   // copy or reference to the original point
    int index = -1;
    su2double distance = 0.0;
};

class CDatapointcloud {
public:
    // --- Data ---
    std::vector<CDatapoint> points;

    // --- Constructors ---
    CDatapointcloud() = default;

    // Construct from vector of points
    explicit CDatapointcloud(const std::vector<CDatapoint>& points) : points(points) {}

    // Construct from matrix of doubles, which are casted into points (coordinates only)
    explicit CDatapointcloud(const Mat& coords) {
        points.reserve(coords.size());
        for (const auto& c : coords)
            points.emplace_back(c);  // uses Point(const Vec&)
    }

    // Construct from coordinates + values
    CDatapointcloud(const Mat& coords, const Mat& vals) {
        if (coords.size() != vals.size()) {
            throw std::runtime_error("Pointcloud: coords and vals size mismatch");
        }
        points.reserve(coords.size());
        for (std::size_t i = 0; i < coords.size(); ++i)
            points.emplace_back(coords[i], vals[i]);
    }

    // --- Accessors ---
    std::size_t size() const { return points.size(); }
    bool empty() const { return points.empty(); }

    // --- Methods ---
    Neighbor nearest(const CDatapoint& query) const;
    std::vector<Neighbor> kNearest(const CDatapoint& query, int k) const;
};
