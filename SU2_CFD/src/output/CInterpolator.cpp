#include "../../include/output/CInterpolator.hpp"
#include <cmath>
#include <limits>
#include <stdexcept>


CDatapoint interpolate_IDW(const CDatapointcloud& pc,
                      const CDatapoint& query_point,
                      int k,
                      su2double power)
{
if (pc.size() == 0) {
    throw std::runtime_error("interpolate_IDW: Pointcloud is empty");
}

    // Find k nearest neighbors
    std::vector<Neighbor> neighbors = pc.kNearest(query_point, k);

    // Compute IDW interpolation
    Vec interpolated_values;
    interpolated_values.resize(pc.points[0].values.size(), 0.0);
    su2double weight_sum = 0.0;

    for (const auto& neighbor : neighbors) {
        su2double dist = neighbor.distance;
        su2double weight;
        if (dist == 0.0) {
            // If the query point coincides with a data point, return its values directly
            return neighbor.point;
        } else {
            weight = 1.0 / std::pow(dist, power);
        }

        for (std::size_t i = 0; i < interpolated_values.size(); ++i) {
            interpolated_values[i] += neighbor.point.values[i] * weight;
        }
        weight_sum += weight;
    }

    // Normalize by total weight
    for (auto& val : interpolated_values) {
        val /= weight_sum;
    }

    return CDatapoint(query_point.coords, interpolated_values);
}

// Single query
CDatapoint Interpolator::interpolate(const CDatapointcloud& pc,
                                 const CDatapoint& query_point,
                                 const std::string& method_in) const
{
    std::string method_to_use = method_in.empty() ? method : method_in;

    if (method_to_use == "IDW") {
        return interpolate_IDW(pc, query_point, k, power);
    } else {
        throw std::runtime_error("Interpolator::interpolate: unknown method '" + method_to_use + "'");
    }
}

// Many queries
CDatapointcloud Interpolator::interpolate(const CDatapointcloud& pc,
                                        const CDatapointcloud& queries,
                                        const std::string& method_in) const
{
    std::vector<CDatapoint> out;
    out.reserve(queries.points.size());
    for (const auto& q : queries.points)
        out.push_back(interpolate(pc, q, method_in));
    return CDatapointcloud(out);
}