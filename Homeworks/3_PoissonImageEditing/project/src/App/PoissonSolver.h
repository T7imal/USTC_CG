#pragma once
#include <opencv2/core/core.hpp>
#include <Eigen/Dense>
#include <Eigen/Sparse>

struct Edge {
    double x;
    int ymax;
    double dx;
};

class PoissonSolver {
public:
    PoissonSolver(std::vector<std::list<Edge>> edge_table, int ymin, cv::Mat source_image, cv::Mat target_image, int pos_x, int pos_y);
    PoissonSolver(Eigen::Vector2i rect_start, Eigen::Vector2i rect_end, cv::Mat source_image, cv::Mat target_image, int pos_x, int pos_y);
    ~PoissonSolver();

    cv::Mat solve();

private:
    bool is_polygon_ = false;
    std::vector<std::list<Edge>> edge_table_;
    bool is_rect_ = false;
    Eigen::Vector2i rect_start_;
    Eigen::Vector2i rect_end_;
    cv::Mat source_image_;
    cv::Mat target_image_;
    // std::vector<bool> inside_mask_;
    Eigen::MatrixXi inside_mask_;
    Eigen::MatrixXi index;
    int num_inside_pixel_ = 0;
    Eigen::Vector2i pos;
};



