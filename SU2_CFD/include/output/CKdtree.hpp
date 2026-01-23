#pragma once

#include <vector>
#include <algorithm>
#include <queue>
#include <limits>
#include <stdexcept>
#include <utility>

#include "CMeshPoint.hpp"
#include "CMeshPointCloud.hpp"

class CKDTree {
private:
    struct KDNode {
        CMeshPoint point;
        KDNode *left = nullptr;
        KDNode *right = nullptr;

        explicit KDNode(CMeshPoint p) : point(std::move(p)) {}
    };

    KDNode* root = nullptr;
    std::size_t nDim = 0;   // runtime dimension

    void destroy(KDNode* node) {
        if (!node) return;
        destroy(node->left);
        destroy(node->right);
        delete node;
    }

    // squared distance using runtime dimension
    su2double dist2(const CMeshPoint& a, const CMeshPoint& b) const {
        su2double sum = 0.0;
        for (std::size_t i = 0; i < nDim; ++i) {
            su2double d = a.coord[i] - b.coord[i];
            sum += d * d;
        }
        return sum;
    }

    KDNode* buildTree(std::vector<CMeshPoint>& pts, std::size_t depth) {
        if (pts.empty()) return nullptr;

        const std::size_t axis = (nDim == 0) ? 0 : (depth % nDim);

        std::sort(pts.begin(), pts.end(),
                  [axis](const CMeshPoint& a, const CMeshPoint& b) {
                      return a.coord[axis] < b.coord[axis];
                  });

        const std::size_t median = pts.size() / 2;

        // Move median point into node
        KDNode* node = new KDNode(std::move(pts[median]));

        std::vector<CMeshPoint> left(pts.begin(), pts.begin() + median);
        std::vector<CMeshPoint> right(pts.begin() + median + 1, pts.end());

        node->left = buildTree(left, depth + 1);
        node->right = buildTree(right, depth + 1);
        return node;
    }

    struct Candidate {
        su2double d2;
        const CMeshPoint* point;  // pointer to avoid copying values during search
    };
    struct WorseFirst {
        bool operator()(const Candidate& a, const Candidate& b) const {
            return a.d2 < b.d2; // max-heap by distance
        }
    };

    void knnSearch(KDNode* node,
                   const CMeshPoint& query,
                   std::size_t k,
                   std::priority_queue<Candidate, std::vector<Candidate>, WorseFirst>& heap,
                   std::size_t depth) const
    {
        if (!node) return;

        const su2double d2_here = dist2(query, node->point);

        if (heap.size() < k) {
            heap.push(Candidate{d2_here, &node->point});
        } else if (d2_here < heap.top().d2) {
            heap.pop();
            heap.push(Candidate{d2_here, &node->point});
        }

        const std::size_t axis = (nDim == 0) ? 0 : (depth % nDim);
        const su2double diff = query.coord[axis] - node->point.coord[axis];

        KDNode* nearChild = (diff < 0) ? node->left : node->right;
        KDNode* farChild  = (diff < 0) ? node->right : node->left;

        knnSearch(nearChild, query, k, heap, depth + 1);

        const su2double diff2 = diff * diff;
        const su2double worstBestD2 = (heap.size() < k)
            ? std::numeric_limits<su2double>::infinity()
            : heap.top().d2;

        if (diff2 <= worstBestD2) {
            knnSearch(farChild, query, k, heap, depth + 1);
        }
    }

public:
    CKDTree() = default;

    explicit CKDTree(std::size_t dim) : nDim(dim) {}

    ~CKDTree() { destroy(root); }

    std::size_t dim() const noexcept { return nDim; }

    // Build from points
    void build(std::vector<CMeshPoint> points) {
        destroy(root);
        root = nullptr;

        if (points.empty()) { nDim = 0; return; }

        // Set dimension from first point
        nDim = points[0].coord.size();
        if (nDim == 0) throw std::runtime_error("CKDTree::build: point has zero dimension");

        // Validate dimensions
        for (const auto& p : points) {
            if (p.coord.size() != nDim) {
                throw std::runtime_error("CKDTree::build: coordinate dimension mismatch");
            }
        }

        root = buildTree(points, 0);
    }

    // Build from cloud
    void build(const CMeshPointCloud& cloud) {
        build(cloud.points);
    }

    // Single query
    std::vector<CMeshPoint> kNearest(const CMeshPoint& query, std::size_t k) const {
        std::vector<CMeshPoint> result;
        if (!root || k == 0) return result;

        if (query.coord.size() != nDim) {
            throw std::runtime_error("CKDTree::kNearest: query dimension mismatch");
        }

        std::priority_queue<Candidate, std::vector<Candidate>, WorseFirst> heap;
        knnSearch(root, query, k, heap, 0);

        result.reserve(heap.size());
        while (!heap.empty()) {
            result.push_back(*heap.top().point); // copy only k results
            heap.pop();
        }
        std::reverse(result.begin(), result.end());
        return result;
    }

    // Batch queries from a cloud
    std::vector<std::vector<CMeshPoint>> kNearest(const CMeshPointCloud& queries, std::size_t k) const {
        std::vector<std::vector<CMeshPoint>> all;
        all.reserve(queries.points.size());

        for (const auto& q : queries.points) {
            all.push_back(kNearest(q, k));
        }
        return all;
    }

    struct NeighborHit {
        const CMeshPoint* point;  // pointer to stored point (do NOT MPI-send this)
        su2double d2;             // squared distance to query
    };

    std::vector<NeighborHit>kNearestWithDist2(const CMeshPoint& query, std::size_t k) const {
        std::vector<NeighborHit> result;
        if (!root || k == 0) return result;

        if (query.coord.size() != nDim) {
            throw std::runtime_error("CKDTree::kNearestWithDist2: query dimension mismatch");
        }

        std::priority_queue<Candidate, std::vector<Candidate>, WorseFirst> heap;
        knnSearch(root, query, k, heap, 0);

        result.reserve(heap.size());
        while (!heap.empty()) {
            const Candidate& c = heap.top();
            result.push_back(NeighborHit{ c.point, c.d2 });
            heap.pop();
        }

        // heap pops worst-first → reverse to get closest-first
        std::reverse(result.begin(), result.end());
        return result;
    }

    std::vector<std::vector<NeighborHit>>kNearestWithDist2(const CMeshPointCloud& queries, std::size_t k) const {
        std::vector<std::vector<NeighborHit>> all;
        all.reserve(queries.points.size());

        for (const auto& q : queries.points) {
            all.push_back(kNearestWithDist2(q, k));
        }
        return all;
    }

};
