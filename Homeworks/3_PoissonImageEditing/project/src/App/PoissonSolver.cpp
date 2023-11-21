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
        for (int i = 0; i < source_image_.rows; ++i) {
            for (int j = 0; j < source_image_.cols; ++j) {
                if (inside_mask_(i, j)) {
                    int m = index(i, j);
                    tripletList.push_back(Eigen::Triplet<double>(m, m, 4)); //不考虑图片边界
                    if (i > rect_start_.y() && inside_mask_(i - 1, j))
                        tripletList.push_back(Eigen::Triplet<double>(m, index(i - 1, j), -1));
                    if (i < rect_end_.y() && inside_mask_(i + 1, j))
                        tripletList.push_back(Eigen::Triplet<double>(m, index(i + 1, j), -1));
                    if (j > rect_start_.x() && inside_mask_(i, j - 1))
                        tripletList.push_back(Eigen::Triplet<double>(m, index(i, i - 1), -1));
                    if (j < rect_end_.x() && inside_mask_(i, j + 1))
                        tripletList.push_back(Eigen::Triplet<double>(m, index(i, j + 1), -1));
                }
            }
        }
        Eigen::SparseMatrix<double> A(num_inside_pixel_, num_inside_pixel_);
        A.setFromTriplets(tripletList.begin(), tripletList.end());
        A.makeCompressed();
        // 构造不同颜色的b
        Eigen::VectorXd b_red(num_inside_pixel_), b_green(num_inside_pixel_), b_blue(num_inside_pixel_);
        b_red.setZero();
        b_green.setZero();
        b_blue.setZero();
        for (int i = 0; i < source_image_.rows; ++i) {
            for (int j = 0; j < source_image_.cols; ++j) {
                if (inside_mask_(i, j)) {
                    int m = index(i, j);
                    b_red[m] = 4 * source_image_.at<cv::Vec3b>(i, j)[0]
                        - source_image_.at<cv::Vec3b>(i + 1, j)[0]
                        - source_image_.at<cv::Vec3b>(i - 1, j)[0]
                        - source_image_.at<cv::Vec3b>(i, j + 1)[0]
                        - source_image_.at<cv::Vec3b>(i, j - 1)[0]; //不考虑图片边界
                    b_green[m] = 4 * source_image_.at<cv::Vec3b>(i, j)[1]
                        - source_image_.at<cv::Vec3b>(i + 1, j)[1]
                        - source_image_.at<cv::Vec3b>(i - 1, j)[1]
                        - source_image_.at<cv::Vec3b>(i, j + 1)[1]
                        - source_image_.at<cv::Vec3b>(i, j - 1)[1]; //不考虑图片边界
                    b_blue[m] = 4 * source_image_.at<cv::Vec3b>(i, j)[2]
                        - source_image_.at<cv::Vec3b>(i + 1, j)[2]
                        - source_image_.at<cv::Vec3b>(i - 1, j)[2]
                        - source_image_.at<cv::Vec3b>(i, j + 1)[2]
                        - source_image_.at<cv::Vec3b>(i, j - 1)[2]; //不考虑图片边界
                    int row = i - rect_start_.y() + pos.y(), col = j - rect_start_.x() + pos.x();   //在目标图像中的坐标
                    if (i == rect_start_.y() || (!inside_mask_(i - 1, j) && i > rect_start_.y())) {
                        b_red[m] += target_image_.at<cv::Vec3b>(row - 1, col)[0];
                        b_green[m] += target_image_.at<cv::Vec3b>(row - 1, col)[1];
                        b_blue[m] += target_image_.at<cv::Vec3b>(row - 1, col)[2];
                    }
                    if (i == rect_end_.y() || (!inside_mask_(i + 1, j) && i < rect_end_.y())) {
                        b_red[m] += target_image_.at<cv::Vec3b>(row + 1, col)[0];
                        b_green[m] += target_image_.at<cv::Vec3b>(row + 1, col)[1];
                        b_blue[m] += target_image_.at<cv::Vec3b>(row + 1, col)[2];
                    }
                    if (j == rect_start_.x() || (!inside_mask_(i, j - 1) && j > rect_start_.x())) {
                        b_red[m] += target_image_.at<cv::Vec3b>(row, col - 1)[0];
                        b_green[m] += target_image_.at<cv::Vec3b>(row, col - 1)[1];
                        b_blue[m] += target_image_.at<cv::Vec3b>(row, col - 1)[2];
                    }
                    if (j == rect_end_.x() || (!inside_mask_(i, j + 1) && j < rect_end_.x())) {
                        b_red[m] += target_image_.at<cv::Vec3b>(row, col + 1)[0];
                        b_green[m] += target_image_.at<cv::Vec3b>(row, col + 1)[1];
                        b_blue[m] += target_image_.at<cv::Vec3b>(row, col + 1)[2];
                    }
                }
            }
        }

        // 解方程
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
        solver.compute(A);
        // Eigen::VectorXd x_red = solver.solve(b_red);
        // Eigen::VectorXd x_green = solver.solve(b_green);
        // Eigen::VectorXd x_blue = solver.solve(b_blue);
        Eigen::VectorXd x_red(num_inside_pixel_);
        Eigen::VectorXd x_green(num_inside_pixel_);
        Eigen::VectorXd x_blue(num_inside_pixel_);
        x_red = solver.solve(b_red);
        x_green = solver.solve(b_green);
        x_blue = solver.solve(b_blue);
        // 生成结果
        cv::Mat result = target_image_.clone();
        for (int i = 0; i < source_image_.rows; ++i) {
            for (int j = 0; j < source_image_.cols; ++j) {
                if (inside_mask_(i, j)) {
                    int m = index(i, j);
                    int row = i - rect_start_.y() + pos.y(), col = j - rect_start_.x() + pos.x();   //在目标图像中的坐标
                    int red = x_red[m], green = x_green[m], blue = x_blue[m];
                    if (row >= 0 && row < target_image_.rows && col >= 0 && col < target_image_.cols) {
                        result.at<cv::Vec3b>(row, col)[0] = red > 255 ? 255 : (red < 0 ? 0 : red);
                        result.at<cv::Vec3b>(row, col)[1] = green > 255 ? 255 : (green < 0 ? 0 : green);
                        result.at<cv::Vec3b>(row, col)[2] = blue > 255 ? 255 : (blue < 0 ? 0 : blue);
                    }
                }
            }
        }
        return result;
    }

    if (is_polygon_) {

    }
}