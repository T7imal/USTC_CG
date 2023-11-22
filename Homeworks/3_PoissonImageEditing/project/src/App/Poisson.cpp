#include "Poisson.h"
#include <ctime>
#include <stdexcept>

Poisson::Poisson() {
	pixels_num_ = 0;
	index_matrix_.resize(0, 0);
	sparse_matrix_.resize(0, 0);
}

Poisson::~Poisson() {

}

void Poisson::set_insidemask(Eigen::MatrixXi inside_mask) {
	inside_mask_ = inside_mask;
}

void Poisson::PoissonInit(cv::Mat source_img) {
	width_ = inside_mask_.rows();
	height_ = inside_mask_.cols();
	paste_point_.setX(0);
	paste_point_.setY(0);
	index_matrix_.resize(width_, height_);
	index_matrix_.setZero();
	for (int i = 0; i < width_; i++) {
		for (int j = 0; j < height_; j++) {
			if (inside_mask_(i, j) == 1) {
				index_matrix_(i, j) = pixels_num_;
				pixels_num_++;
			}
		}
	}

	sparse_matrix_.resize(pixels_num_, pixels_num_);
	sparse_matrix_.setZero();

	QVector<Eigen::Triplet<float>> tripletList;
	for (int i = 0; i < width_; i++) {
		for (int j = 0; j < height_; j++) {
			if (inside_mask_(i, j) == 1) {
				int index = index_matrix_(i, j);
				tripletList.push_back(Eigen::Triplet<float>(index, index, 4));
				if (i > 0 && inside_mask_(i - 1, j) == 1) {
					tripletList.push_back(Eigen::Triplet<float>(index, index_matrix_(i - 1, j), -1));
				}
				if (j > 0 && inside_mask_(i, j - 1) == 1) {
					tripletList.push_back(Eigen::Triplet<float>(index, index_matrix_(i, j - 1), -1));
				}
				if (i < width_ - 1 && inside_mask_(i + 1, j) == 1) {
					tripletList.push_back(Eigen::Triplet<float>(index, index_matrix_(i + 1, j), -1));
				}
				if (j < height_ - 1 && inside_mask_(i, j + 1) == 1) {
					tripletList.push_back(Eigen::Triplet<float>(index, index_matrix_(i, j + 1), -1));
				}
			}
		}
	}
	sparse_matrix_.setFromTriplets(tripletList.begin(), tripletList.end());
	sparse_matrix_.makeCompressed();

	solver.compute(sparse_matrix_);
}

void Poisson::GetPoisson(QPoint paste_point, QPoint source_point, cv::Mat& paste_img_, cv::Mat& source_img_) {
	paste_point_ = paste_point;
	source_point_ = source_point;

	div_red_.resize(pixels_num_);
	div_green_.resize(pixels_num_);
	div_blue_.resize(pixels_num_);
	div_red_.setZero();
	div_green_.setZero();
	div_blue_.setZero();

	for (int i = 0; i < width_; i++) {
		for (int j = 0; j < height_; j++) {
			if (inside_mask_(i, j) == 1) {
				int index = index_matrix_(i, j);
				int x = source_point_.y() + i;
				int y = source_point_.x() + j;

				div_red_[index] = 4 * source_img_.at<cv::Vec3b>(x, y)[0]
					- source_img_.at<cv::Vec3b>(x + 1, y)[0]
					- source_img_.at<cv::Vec3b>(x - 1, y)[0]
					- source_img_.at<cv::Vec3b>(x, y + 1)[0]
					- source_img_.at<cv::Vec3b>(x, y - 1)[0]; //不考虑图片边界
				div_green_[index] = 4 * source_img_.at<cv::Vec3b>(x, y)[1]
					- source_img_.at<cv::Vec3b>(x + 1, y)[1]
					- source_img_.at<cv::Vec3b>(x - 1, y)[1]
					- source_img_.at<cv::Vec3b>(x, y + 1)[1]
					- source_img_.at<cv::Vec3b>(x, y - 1)[1]; //不考虑图片边界
				div_blue_[index] = 4 * source_img_.at<cv::Vec3b>(x, y)[2]
					- source_img_.at<cv::Vec3b>(x + 1, y)[2]
					- source_img_.at<cv::Vec3b>(x - 1, y)[2]
					- source_img_.at<cv::Vec3b>(x, y + 1)[2]
					- source_img_.at<cv::Vec3b>(x, y - 1)[2]; //不考虑图片边界

				if (i == 0 || (i > 0 && inside_mask_(i - 1, j) == 0)) {
					div_red_[index] += paste_img_.at<cv::Vec3b>(i + paste_point_.y() - 1, j + paste_point_.x())[0];
					div_green_[index] += paste_img_.at<cv::Vec3b>(i + paste_point_.y() - 1, j + paste_point_.x())[1];
					div_blue_[index] += paste_img_.at<cv::Vec3b>(i + paste_point_.y() - 1, j + paste_point_.x())[2];
				}
				if (i == width_ - 1 || (i < width_ - 1 && inside_mask_(i + 1, j) == 0)) {
					div_red_[index] += paste_img_.at<cv::Vec3b>(i + paste_point_.y() + 1, j + paste_point_.x())[0];
					div_green_[index] += paste_img_.at<cv::Vec3b>(i + paste_point_.y() + 1, j + paste_point_.x())[1];
					div_blue_[index] += paste_img_.at<cv::Vec3b>(i + paste_point_.y() + 1, j + paste_point_.x())[2];
				}
				if (j == 0 || (j > 0 && inside_mask_(i, j - 1) == 0)) {
					div_red_[index] += paste_img_.at<cv::Vec3b>(i + paste_point_.y(), j + paste_point_.x() - 1)[0];
					div_green_[index] += paste_img_.at<cv::Vec3b>(i + paste_point_.y(), j + paste_point_.x() - 1)[1];
					div_blue_[index] += paste_img_.at<cv::Vec3b>(i + paste_point_.y(), j + paste_point_.x() - 1)[2];
				}
				if (j == height_ - 1 || (j < height_ - 1 && inside_mask_(i, j + 1) == 0)) {
					div_red_[index] += paste_img_.at<cv::Vec3b>(i + paste_point_.y(), j + paste_point_.x() + 1)[0];
					div_green_[index] += paste_img_.at<cv::Vec3b>(i + paste_point_.y(), j + paste_point_.x() + 1)[1];
					div_blue_[index] += paste_img_.at<cv::Vec3b>(i + paste_point_.y(), j + paste_point_.x() + 1)[2];
				}
			}
		}
	}

	Eigen::VectorXf vec_red(pixels_num_);
	Eigen::VectorXf vec_green(pixels_num_);
	Eigen::VectorXf vec_blue(pixels_num_);

	vec_red = solver.solve(div_red_);
	vec_green = solver.solve(div_green_);
	vec_blue = solver.solve(div_blue_);

	for (int i = 0; i < width_; i++) {
		for (int j = 0; j < height_; j++) {
			if (inside_mask_(i, j) == 1) {
				int index = index_matrix_(i, j);
				int red = vec_red(index);
				int green = vec_green(index);
				int blue = vec_blue(index);
				int x = i + paste_point_.y(), y = j + paste_point_.x();
				paste_img_.at<cv::Vec3b>(x, y)[0] = red > 255 ? 255 : (red < 0 ? 0 : red);
				paste_img_.at<cv::Vec3b>(x, y)[1] = green > 255 ? 255 : (green < 0 ? 0 : green);
				paste_img_.at<cv::Vec3b>(x, y)[2] = blue > 255 ? 255 : (blue < 0 ? 0 : blue);
			}
		}
	}
}

double Poisson::VecLength(cv::Vec3i vec) {
	return sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
}

void Poisson::MixingPoisson(QPoint paste_point, QPoint source_point, cv::Mat& paste_img_, cv::Mat& source_img_) {
	paste_point_ = paste_point;
	source_point_ = source_point;

	div_red_.resize(pixels_num_);
	div_green_.resize(pixels_num_);
	div_blue_.resize(pixels_num_);
	div_red_.setZero();
	div_green_.setZero();
	div_blue_.setZero();

	for (int i = 0; i < width_; i++) {
		for (int j = 0; j < height_; j++) {
			if (inside_mask_(i, j) == 1) {
				int index = index_matrix_(i, j);
				int x = source_point_.y() + i;
				int y = source_point_.x() + j;
				int x_ = paste_point_.y() + i;
				int y_ = paste_point_.x() + j;
				cv::Vec3i vec, temp_vec, temp_vec_paste;

				temp_vec = source_img_.at<cv::Vec3b>(x, y);
				temp_vec -= source_img_.at<cv::Vec3b>(x + 1, y);
				temp_vec_paste = paste_img_.at<cv::Vec3b>(x_, y_);
				temp_vec_paste -= paste_img_.at<cv::Vec3b>(x_ + 1, y_);
				vec = VecLength(temp_vec) > VecLength(temp_vec_paste) ? temp_vec : temp_vec_paste;

				temp_vec = source_img_.at<cv::Vec3b>(x, y);
				temp_vec -= source_img_.at<cv::Vec3b>(x - 1, y);
				temp_vec_paste = paste_img_.at<cv::Vec3b>(x_, y_);
				temp_vec_paste -= paste_img_.at<cv::Vec3b>(x_ - 1, y_);
				vec += VecLength(temp_vec) > VecLength(temp_vec_paste) ? temp_vec : temp_vec_paste;

				temp_vec = source_img_.at<cv::Vec3b>(x, y);
				temp_vec -= source_img_.at<cv::Vec3b>(x, y + 1);
				temp_vec_paste = paste_img_.at<cv::Vec3b>(x_, y_);
				temp_vec_paste -= paste_img_.at<cv::Vec3b>(x_, y_ + 1);
				vec += VecLength(temp_vec) > VecLength(temp_vec_paste) ? temp_vec : temp_vec_paste;

				temp_vec = source_img_.at<cv::Vec3b>(x, y);
				temp_vec -= source_img_.at<cv::Vec3b>(x, y - 1);
				temp_vec_paste = paste_img_.at<cv::Vec3b>(x_, y_);
				temp_vec_paste -= paste_img_.at<cv::Vec3b>(x_, y_ - 1);
				vec += VecLength(temp_vec) > VecLength(temp_vec_paste) ? temp_vec : temp_vec_paste;

				div_red_(index_matrix_(i, j)) = vec[0];
				div_green_(index_matrix_(i, j)) = vec[1];
				div_blue_(index_matrix_(i, j)) = vec[2];


				if (i == 0 || (i > 0 && inside_mask_(i - 1, j) == 0)) {
					div_red_[index] += paste_img_.at<cv::Vec3b>(i + paste_point_.y() - 1, j + paste_point_.x())[0];
					div_green_[index] += paste_img_.at<cv::Vec3b>(i + paste_point_.y() - 1, j + paste_point_.x())[1];
					div_blue_[index] += paste_img_.at<cv::Vec3b>(i + paste_point_.y() - 1, j + paste_point_.x())[2];
				}
				if (i == width_ - 1 || (i < width_ - 1 && inside_mask_(i + 1, j) == 0)) {
					div_red_[index] += paste_img_.at<cv::Vec3b>(i + paste_point_.y() + 1, j + paste_point_.x())[0];
					div_green_[index] += paste_img_.at<cv::Vec3b>(i + paste_point_.y() + 1, j + paste_point_.x())[1];
					div_blue_[index] += paste_img_.at<cv::Vec3b>(i + paste_point_.y() + 1, j + paste_point_.x())[2];
				}
				if (j == 0 || (j > 0 && inside_mask_(i, j - 1) == 0)) {
					div_red_[index] += paste_img_.at<cv::Vec3b>(i + paste_point_.y(), j + paste_point_.x() - 1)[0];
					div_green_[index] += paste_img_.at<cv::Vec3b>(i + paste_point_.y(), j + paste_point_.x() - 1)[1];
					div_blue_[index] += paste_img_.at<cv::Vec3b>(i + paste_point_.y(), j + paste_point_.x() - 1)[2];
				}
				if (j == height_ - 1 || (j < height_ - 1 && inside_mask_(i, j + 1) == 0)) {
					div_red_[index] += paste_img_.at<cv::Vec3b>(i + paste_point_.y(), j + paste_point_.x() + 1)[0];
					div_green_[index] += paste_img_.at<cv::Vec3b>(i + paste_point_.y(), j + paste_point_.x() + 1)[1];
					div_blue_[index] += paste_img_.at<cv::Vec3b>(i + paste_point_.y(), j + paste_point_.x() + 1)[2];
				}
			}
		}
	}

	Eigen::VectorXf vec_red(pixels_num_);
	Eigen::VectorXf vec_green(pixels_num_);
	Eigen::VectorXf vec_blue(pixels_num_);

	vec_red = solver.solve(div_red_);
	vec_green = solver.solve(div_green_);
	vec_blue = solver.solve(div_blue_);

	for (int i = 0; i < width_; i++) {
		for (int j = 0; j < height_; j++) {
			if (inside_mask_(i, j) == 1) {
				int index = index_matrix_(i, j);
				int red = vec_red(index);
				int green = vec_green(index);
				int blue = vec_blue(index);
				int x = i + paste_point_.y(), y = j + paste_point_.x();
				paste_img_.at<cv::Vec3b>(x, y)[0] = red > 255 ? 255 : (red < 0 ? 0 : red);
				paste_img_.at<cv::Vec3b>(x, y)[1] = green > 255 ? 255 : (green < 0 ? 0 : green);
				paste_img_.at<cv::Vec3b>(x, y)[2] = blue > 255 ? 255 : (blue < 0 ? 0 : blue);
			}
		}
	}
}

void Poisson::Paste(QPoint paste_point, QPoint source_point, cv::Mat& paste_img_, cv::Mat& source_img_) {
	paste_point_ = paste_point;
	source_point_ = source_point;
	for (int j = 0; j < width_; j++) {
		for (int k = 0; k < height_; k++) {
			if (inside_mask_(j, k) == 1) {
				paste_img_.at<cv::Vec3b>(j + paste_point_.y(), k + paste_point_.x())[0] = source_img_.at<cv::Vec3b>(j + source_point_.y(), k + source_point_.x())[0];
				paste_img_.at<cv::Vec3b>(j + paste_point_.y(), k + paste_point_.x())[1] = source_img_.at<cv::Vec3b>(j + source_point_.y(), k + source_point_.x())[1];
				paste_img_.at<cv::Vec3b>(j + paste_point_.y(), k + paste_point_.x())[2] = source_img_.at<cv::Vec3b>(j + source_point_.y(), k + source_point_.x())[2];
			}
		}
	}
}