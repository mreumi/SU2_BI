#pragma once

#include <vector>
#include <cstddef>
#include <stdexcept>

#include "CMeshPoint.hpp"

// Type aliases (same as before)
using Vec = std::vector<su2double>;
using Mat = std::vector<Vec>;

/**
 * CMeshPointCloud<Dim>
 *
 * - Owns a collection of CMeshPoint<Dim>
 */
class CMeshPointCloud {
public:
    // --- Data ---
    std::vector<CMeshPoint> points;

    // --- Constructors ---
    CMeshPointCloud() = default;

    explicit CMeshPointCloud(const std::vector<CMeshPoint>& pts)
        : points(pts) {}

    /**
     * Construct from coordinate matrix only.
     * Values are initialized to zero.
     *
     * coords[i].size() must be exactly Dim.
     */
    explicit CMeshPointCloud(const Mat& coords) {
        points.reserve(coords.size());

        for (const auto& c : coords) {
            if (c.size() != coords[0].size()) {
                throw std::runtime_error(
                    "CMeshPointCloud: coordinate dimension mismatch"
                );
            }

            std::vector<su2double> a(c.size());
            for (std::size_t d = 0; d < c.size(); ++d) {
                a[d] = c[d];
            }

            points.emplace_back(a, std::vector<su2double>{});
        }
    }

    CMeshPointCloud(const Mat& coords, const Mat& vals) {
        if (coords.size() != vals.size()) {
            throw std::runtime_error(
                "CMeshPointCloud: coords and vals size mismatch"
            );
        }

        points.reserve(coords.size());

        for (std::size_t i = 0; i < coords.size(); ++i) {
            if (coords[i].size() != coords[0].size()) {
                throw std::runtime_error(
                    "CMeshPointCloud: coordinate dimension mismatch"
                );
            }

            std::vector<su2double> a(coords[i].size());
            for (std::size_t d = 0; d < coords[i].size(); ++d) {
                a[d] = coords[i][d];
            }

            // Store the FULL value vector
            points.emplace_back(a, vals[i]);
        }
    }

    // --- Accessors ---

    std::size_t size() const noexcept { return points.size(); }
    bool empty() const noexcept { return points.empty(); }

    const std::vector<CMeshPoint>& get_points() const noexcept {
        return points;
    }
};
