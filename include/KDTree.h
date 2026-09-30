#ifndef PAM_ST_KDTREE_H
#define PAM_ST_KDTREE_H

#include <vector>

// KD-tree over weighted points in R^dim for fixed-radius neighbour counting.
class KDTreeND {
public:
    // points holds weights.size() rows of dim coordinates; both must outlive the tree.
    KDTreeND(const std::vector<double>& points, int dim, const std::vector<int>& weights);

    // counts[i] = total weight of points j with squared distance(i, j) <= radius2,
    // including i itself. Each pair is tested once and credited to both points,
    // which the symmetry of the distance test makes exact.
    std::vector<int> weighted_neighbor_counts(double radius2, int threads) const;

    // Total weight of points within radius2 of an arbitrary query row, which
    // need not be one of the tree's own points.
    int weighted_count_within(const double* query, double radius2) const;

private:
    struct Node {
        int lo = 0;
        int hi = 0;
        int left = -1;
        int right = -1;
    };

    // Counts a query contributes to itself and to the points it is paired with.
    struct PairCounts {
        std::vector<int> own;     // per point, from pairs found by that point's own query
        std::vector<int> ranges;  // difference array over whole subtrees credited at once
    };

    static constexpr int kLeafSize = 16;

    const std::vector<double>& points_;
    const std::vector<int>& weights_;
    int n_ = 0;
    int dim_ = 0;
    std::vector<int> indices_;        // input index of each ordered position
    std::vector<Node> nodes_;
    std::vector<double> bbox_min_;
    std::vector<double> bbox_max_;
    std::vector<double> ordered_points_;
    std::vector<int> ordered_weights_;
    std::vector<int> weight_prefix_;  // weight_prefix_[i] = total weight of positions < i

    const double* row(int point) const;
    const double* ordered_row(int ordered_pos) const;
    int build(int lo, int hi);
    bool within_radius2(int ordered_pos, const double* query, double radius2) const;
    bool bbox_outside_radius2(int node_idx, const double* query, double radius2) const;
    bool bbox_inside_radius2(int node_idx, const double* query, double radius2) const;
    void count_pairs_above(int node_idx, int ordered_pos, const double* query,
                           double radius2, PairCounts& counts) const;
    int count_within(int node_idx, const double* query, double radius2) const;
};

#endif  // PAM_ST_KDTREE_H
