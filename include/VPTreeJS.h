#ifndef PAM_ST_VPTREEJS_H
#define PAM_ST_VPTREEJS_H

#include <vector>

// Vantage-point tree over weighted probability vectors, using the square root
// of the Jensen-Shannon divergence as the metric.
class VPTreeJS {
public:
    // points holds weights.size() rows of dim probabilities with their Shannon
    // entropies; all three vectors must outlive the tree.
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
    // Vantage points copied per node, so a query walks contiguous memory.
    std::vector<double> node_points_;
    std::vector<double> node_entropies_;

    const double* row(int point) const;
    double distance2(int a, int b) const;
    int select_vantage(int lo, int hi) const;
    int build(int lo, int hi);
    // point is the querying tree point, or -1 for an external query row.
    int count(int node_idx, int point, const double* query, double query_entropy,
              double radius2, double radius) const;
};

#endif  // PAM_ST_VPTREEJS_H
