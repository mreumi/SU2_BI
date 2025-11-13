#include "../../include/output/CDatapointcloud.hpp"

#include <stdexcept>
#include <cmath>
#include <algorithm>  // sort
#include <numeric>    // iota

Neighbor CDatapointcloud::nearest(const CDatapoint& point_in) const {
    if (points.empty())
        throw std::runtime_error("Pointcloud::nearest: no points");

    // Neighbor best is initialized with the first point (point, index, distance)
    Neighbor closest_point{points[0], 0, point_in.distance_to(points[0])};

    // Iterate through all points to find the nearest
    for (std::size_t i = 1; i < points.size(); ++i) {
        su2double d = point_in.distance_to(points[i]);
        if (d < closest_point.distance) {
            closest_point = Neighbor{points[i], static_cast<int>(i), d};
        }
    }
    return closest_point;
}

std::vector<Neighbor> CDatapointcloud::kNearest(const CDatapoint& point_in, int k) const {
    if (points.empty())
        throw std::runtime_error("Pointcloud::kNearest: no points");
    if (k <= 0)
        throw std::runtime_error("Pointcloud::kNearest: k must be positive");

    // Number of points
    const std::size_t n = points.size();
    if (k > static_cast<int>(n))
        k = static_cast<int>(n);

    // Iterate through all points to find the k nearest
    Vec distances;
    std::vector<int> indices;
    distances.reserve(n);
    indices.reserve(n);
    
    for (std::size_t i = 0; i < points.size(); ++i) {
        su2double d = point_in.distance_to(points[i]);
        distances.push_back(d);
        indices.push_back(static_cast<int>(i));
    }

    // Sort distances and indices together
    std::vector<std::size_t> order(indices.size());
    std::iota(order.begin(), order.end(), 0); // fill with 0, 1, ..., n-1

    std::sort(order.begin(), order.end(),
          [&](std::size_t a, std::size_t b) { return distances[a] < distances[b]; });

    // Collect the k nearest neighbors
    std::vector<Neighbor> neighbors;
    for (int i = 0; i < k; ++i) {
        int idx = order[i];
        neighbors.push_back({points[idx], idx, distances[idx]});
    }

    return neighbors;
}

