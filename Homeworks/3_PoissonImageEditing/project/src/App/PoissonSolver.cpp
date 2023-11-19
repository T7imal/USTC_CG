#include "PoissonSolver.h"
#include <iostream>

PoissonSolver::PoissonSolver(std::vector<std::list<Edge>> edge_table, int ymin, cv::Mat source_image, cv::Mat target_image, int pos_x, int pos_y) {
    is_polygon_ = true;
    edge_table_ = edge_table;
    source_image_ = source_image;
    target_image_ = target_image;
    inside_mask_ = Eigen::MatrixXi::Zero(source_image.rows, source_image.cols);
    index = Eigen::MatrixXi::Zero(source_image.rows, source_image.cols);
    pos = Eigen::Vector2i(pos_x, pos_y);
    // 扫描线填充
    // 建立活动边表
    std::list<Edge> active_edge_table;
    auto tmp = edge_table[0].begin();
    while (tmp != edge_table[0].end()) {
        active_edge_table.push_back(*tmp);
        ++tmp;
    }
    std::vector<Eigen::Triplet<int>> tripletList;
    for (int i = 0; i < edge_table.size(); ++i) {
        tmp = active_edge_table.begin();
        while (tmp != active_edge_table.end()) {
            int x1 = tmp->x;
            int x2 = (++tmp)->x;
            for (int j = x1; j <= x2; ++j) {
                inside_mask_(i + ymin, j) = true;
                index(i + ymin, j) = num_inside_pixel_;
                num_inside_pixel_++;
            }
            ++tmp;
        }
        // 删去活动边表中不再相交的边，更新活动边表中的x
        tmp = active_edge_table.begin();
        while (tmp != active_edge_table.end()) {
            if (tmp->ymax == i + 1 + ymin)
                tmp = active_edge_table.erase(tmp);
            else
                ++tmp;
        }
        tmp = active_edge_table.begin();
        while (tmp != active_edge_table.end()) {
            tmp->x += tmp->dx;
            ++tmp;
        }
        // 将新的边加入活动边表
        if (i == edge_table.size() - 1)
            break;
        tmp = edge_table[i + 1].begin();
        while (tmp != edge_table[i + 1].end()) {
            active_edge_table.push_back(*tmp);
            ++tmp;
        }
        active_edge_table.sort([](const Edge& e1, const Edge& e2) {return e1.x < e2.x; });
    }
}

PoissonSolver::PoissonSolver(Eigen::Vector2i rect_start, Eigen::Vector2i rect_end, cv::Mat source_image, cv::Mat target_image, int pos_x, int pos_y) {
    is_rect_ = true;
    rect_start_ = rect_start;
    rect_end_ = rect_end;
    source_image_ = source_image;
    target_image_ = target_image;
    inside_mask_ = Eigen::MatrixXi::Zero(source_image.rows, source_image.cols);
    index = Eigen::MatrixXi::Zero(source_image.rows, source_image.cols);
    pos = Eigen::Vector2i(pos_x, pos_y);
    for (int i = rect_start_.y(); i < rect_end_.y(); ++i) {
        for (int j = rect_start_.x(); j < rect_end_.x(); ++j) {
            inside_mask_(i, j) = true;
            index(i, j) = num_inside_pixel_;
            num_inside_pixel_++;
        }
    }
}

PoissonSolver::~PoissonSolver() {

}

cv::Mat PoissonSolver::solve() {
    if (is_rect_) {
        // 构造系数矩阵
        std::vector<Eigen::Triplet<double>> tripletList;
        for (int i = rect_start_.y(); i < rect_end_.y(); ++i) {
            for (int j = rect_start_.x(); j < rect_end_.x(); ++j) {
                if (inside_mask_(i, j)) {
                    int m = index(i, j);
                    tripletList.push_back(Eigen::Triplet<double>(m, m, 4)); //不考虑图片边界
                    if (inside_mask_(i - 1, j))
                        tripletList.push_back(Eigen::Triplet<double>(m, index(i - 1, j), -1));
                    if (inside_mask_(i + 1, j))
                        tripletList.push_back(Eigen::Triplet<double>(m, index(i + 1, j), -1));
                    if (inside_mask_(i, j - 1))
                        tripletList.push_back(Eigen::Triplet<double>(m, index(i, i - 1), -1));
                    if (inside_mask_(i, j + 1))
                        tripletList.push_back(Eigen::Triplet<double>(m, index(i, j + 1), -1));
                }
            }
        }
        Eigen::SparseMatrix<double> A(num_inside_pixel_, num_inside_pixel_);
        A.setFromTriplets(tripletList.begin(), tripletList.end());
        A.makeCompressed();
        // 构造不同颜色的b
        Eigen::VectorXd b_red(num_inside_pixel_), b_green(num_inside_pixel_), b_blue(num_inside_pixel_);
        for (int i = rect_start_.y(); i < rect_end_.y(); ++i) {
            for (int j = rect_start_.x(); j < rect_end_.x(); ++j) {
                if (inside_mask_(i, j)) {
                    int m = index(i, j);
                    b_red[m] = 4 * target_image_.at<cv::Vec3b>(i + pos.y(), j + pos.x())[2]
                        - target_image_.at<cv::Vec3b>(i + pos.y() - 1, j + pos.x())[2]
                        - target_image_.at<cv::Vec3b>(i + pos.y() + 1, j + pos.x())[2]
                        - target_image_.at<cv::Vec3b>(i + pos.y(), j + pos.x() - 1)[2]
                        - target_image_.at<cv::Vec3b>(i + pos.y(), j + pos.x() + 1)[2]; //不考虑图片边界
                    b_green[m] = 4 * target_image_.at<cv::Vec3b>(i + pos.y(), j + pos.x())[1]
                        - target_image_.at<cv::Vec3b>(i + pos.y() - 1, j + pos.x())[1]
                        - target_image_.at<cv::Vec3b>(i + pos.y() + 1, j + pos.x())[1]
                        - target_image_.at<cv::Vec3b>(i + pos.y(), j + pos.x() - 1)[1]
                        - target_image_.at<cv::Vec3b>(i + pos.y(), j + pos.x() + 1)[1]; //不考虑图片边界
                    b_blue[m] = 4 * target_image_.at<cv::Vec3b>(i + pos.y(), j + pos.x())[0]
                        - target_image_.at<cv::Vec3b>(i + pos.y() - 1, j + pos.x())[0]
                        - target_image_.at<cv::Vec3b>(i + pos.y() + 1, j + pos.x())[0]
                        - target_image_.at<cv::Vec3b>(i + pos.y(), j + pos.x() - 1)[0]
                        - target_image_.at<cv::Vec3b>(i + pos.y(), j + pos.x() + 1)[0]; //不考虑图片边界
                    if (!inside_mask_(i - 1, j)) {
                        b_red[m] += source_image_.at<cv::Vec3b>(i - 1, j)[2];
                        b_green[m] += source_image_.at<cv::Vec3b>(i - 1, j)[1];
                        b_blue[m] += source_image_.at<cv::Vec3b>(i - 1, j)[0];
                    }
                    if (!inside_mask_(i + 1, j)) {
                        b_red[m] += source_image_.at<cv::Vec3b>(i + 1, j)[2];
                        b_green[m] += source_image_.at<cv::Vec3b>(i + 1, j)[1];
                        b_blue[m] += source_image_.at<cv::Vec3b>(i + 1, j)[0];
                    }
                    if (!inside_mask_(i, j - 1)) {
                        b_red[m] += source_image_.at<cv::Vec3b>(i, j - 1)[2];
                        b_green[m] += source_image_.at<cv::Vec3b>(i, j - 1)[1];
                        b_blue[m] += source_image_.at<cv::Vec3b>(i, j - 1)[0];
                    }
                    if (!inside_mask_(i, j + 1)) {
                        b_red[m] += source_image_.at<cv::Vec3b>(i, j + 1)[2];
                        b_green[m] += source_image_.at<cv::Vec3b>(i, j + 1)[1];
                        b_blue[m] += source_image_.at<cv::Vec3b>(i, j + 1)[0];
                    }
                }
            }
        }
        // 解方程
        Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> solver;
        solver.compute(A);
        Eigen::VectorXd x_red = solver.solve(b_red);
        Eigen::VectorXd x_green = solver.solve(b_green);
        Eigen::VectorXd x_blue = solver.solve(b_blue);
        // 生成结果
        cv::Mat result = target_image_.clone();
        for (int i = rect_start_.y(); i < rect_end_.y(); ++i) {
            for (int j = rect_start_.x(); j < rect_end_.x(); ++j) {
                if (inside_mask_(i, j)) {
                    int m = index(i, j);
                    result.at<cv::Vec3b>(i - rect_start_.y() + pos.y(), j - rect_start_.x() + pos.x())[2] = x_red[m];
                    result.at<cv::Vec3b>(i - rect_start_.y() + pos.y(), j - rect_start_.x() + pos.x())[1] = x_green[m];
                    result.at<cv::Vec3b>(i - rect_start_.y() + pos.y(), j - rect_start_.x() + pos.x())[0] = x_blue[m];

                }
            }
        }
        return result;
    }

    if (is_polygon_) {

    }
}