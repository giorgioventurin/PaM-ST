#ifndef PAM_ST_VPTREEJS_H
#define PAM_ST_VPTREEJS_H

#include <vector>

// Vantage-point tree over weighted probability vectors, using the square root
// of the Jensen-Shannon divergence as the metric.
class VPTreeJS {
public:
    // points holds weights.size() rows of dim probabilities with their Shannon
    // entropies, as shannon_entropy computes them; all three vectors must
    // outlive the tree. Query entropies are expected the same way.
    VPTreeJS(const std::vector<double>& points,
             const std::vector<double>& entropies,
             int dim,
             const std::vector<int>& weights);

    // counts[i] = total weight of points j within JS distance `radius` of point i.
    std::vector<int> weighted_neighbor_counts(double radius, int threads) const;

    // Total weight of points within `radius` of an arbitrary probability row
    // and its entropy; the row need not be one of the tree's own points.
    int weighted_count_within(const double* query, double query_entropy, double radius) const;

private:
    struct Node {
        int point = -1;
        int lo = 0;
        int hi = 0;
        int left = -1;
        int right = -1;
        double threshold = 0.0;   // median distance from the vantage point
        double max_radius = 0.0;  // furthest point in the subtree
    };

    const std::vector<double>& points_;
    const std::vector<double>& entropies_;
    const std::vector<int>& weights_;
    int n_ = 0;
    int dim_ = 0;
    std::vector<int> indices_;
    std::vector<Node> nodes_;
    std::vector<int> node_weight_sums_;
    // Mixture-entropy terms each point fixes by itself (see js_row_terms), and
    // the square roots of its probabilities.
    std::vector<double> half_terms_;
    std::vector<double> self_terms_;
    std::vector<double> roots_;
    // Vantage points copied per node, so a query walks contiguous memory.
    std::vector<double> node_points_;
    std::vector<double> node_half_terms_;
    std::vector<double> node_roots_;
    std::vector<double> node_entropies_;
    // Per node, the box spanned by the square roots of every point in its subtree.
    std::vector<double> box_lo_;
    std::vector<double> box_hi_;

    // One query row with its entropy, its own terms and its square roots.
    struct Query {
        const double* row;
        const double* half;
        const double* self;
        const double* roots;
        double entropy;
    };

    const double* row(int point) const;
    const double* half_row(int point) const;
    const double* self_row(int point) const;
    const double* root_row(int point) const;
    double distance2(int a, int b) const;
    int select_vantage(int lo, int hi) const;
    int build(int lo, int hi);
    // point is the querying tree point, or -1 for an external query row.
    int count(int node_idx, int point, const Query& query, double radius2, double radius) const;
};

#endif  // PAM_ST_VPTREEJS_H
