#include "ImageWidget.h"
#include <QImage>
#include <QPainter>
#include <QtWidgets> 
#include <iostream>
#include <vector>
#include <ANN/ANN.h>					// ANN declarations

using std::cout;
using std::endl;

ImageWidget::ImageWidget(void) {
	ptr_image_ = new QImage();
	ptr_image_backup_ = new QImage();
}


ImageWidget::~ImageWidget(void) {
}

void ImageWidget::paintEvent(QPaintEvent* paintevent) {
	QPainter painter;
	painter.begin(this);

	// Draw background
	painter.setBrush(Qt::lightGray);
	QRect back_rect(0, 0, width(), height());
	painter.drawRect(back_rect);

	// Draw image
	QRect rect = QRect((width() - ptr_image_->width()) / 2, (height() - ptr_image_->height()) / 2, ptr_image_->width(), ptr_image_->height());
	painter.drawImage(rect, *ptr_image_);

	if (draw_status_) {
		painter.setPen(Qt::red);
		painter.drawRect(start_point_.x() - 5, start_point_.y() - 5, 10, 10);
		painter.drawRect(end_point_.x() - 5, end_point_.y() - 5, 10, 10);
		painter.drawLine(start_point_, end_point_);
	}

	painter.setPen(Qt::black);
	for (QRect rect : rect_list_) {

		painter.drawRect(rect);
	}
	for (QLine line : line_list_) {
		painter.drawLine(line);
	}

	painter.end();
	update();
}

void ImageWidget::Open() {
	// Open file
	QString fileName = QFileDialog::getOpenFileName(this, tr("Read Image"), ".", tr("Images(*.bmp *.png *.jpg)"));

	// Load file
	if (!fileName.isEmpty()) {
		ptr_image_->load(fileName);
		*(ptr_image_backup_) = *(ptr_image_);
	}

	//ptr_image_->invertPixels(QImage::InvertRgb);
	//*(ptr_image_) = ptr_image_->mirrored(true, true);
	//*(ptr_image_) = ptr_image_->rgbSwapped();
	cout << "image size: " << ptr_image_->width() << ' ' << ptr_image_->height() << endl;
	update();
}

void ImageWidget::Save() {
	SaveAs();
}

void ImageWidget::SaveAs() {
	QString filename = QFileDialog::getSaveFileName(this, tr("Save Image"), ".", tr("Images(*.bmp *.png *.jpg)"));
	if (filename.isNull()) {
		return;
	}

	ptr_image_->save(filename);
}

void ImageWidget::Invert() {
	for (int i = 0; i < ptr_image_->width(); i++) {
		for (int j = 0; j < ptr_image_->height(); j++) {
			QRgb color = ptr_image_->pixel(i, j);
			ptr_image_->setPixel(i, j, qRgb(255 - qRed(color), 255 - qGreen(color), 255 - qBlue(color)));
		}
	}

	// equivalent member function of class QImage
	// ptr_image_->invertPixels(QImage::InvertRgb);
	update();
}

void ImageWidget::Mirror(bool ishorizontal, bool isvertical) {
	QImage image_tmp(*(ptr_image_));
	int width = ptr_image_->width();
	int height = ptr_image_->height();

	if (ishorizontal) {
		if (isvertical) {
			for (int i = 0; i < width; i++) {
				for (int j = 0; j < height; j++) {
					ptr_image_->setPixel(i, j, image_tmp.pixel(width - 1 - i, height - 1 - j));
				}
			}
		}
		else			//仅水平翻转			
		{
			for (int i = 0; i < width; i++) {
				for (int j = 0; j < height; j++) {
					ptr_image_->setPixel(i, j, image_tmp.pixel(width - 1 - i, j));
				}
			}
		}

	}
	else {
		if (isvertical)		//仅垂直翻转
		{
			for (int i = 0; i < width; i++) {
				for (int j = 0; j < height; j++) {
					ptr_image_->setPixel(i, j, image_tmp.pixel(i, height - 1 - j));
				}
			}
		}
	}

	// equivalent member function of class QImage
	//*(ptr_image_) = ptr_image_->mirrored(true, true);
	update();
}

void ImageWidget::TurnGray() {
	for (int i = 0; i < ptr_image_->width(); i++) {
		for (int j = 0; j < ptr_image_->height(); j++) {
			QRgb color = ptr_image_->pixel(i, j);
			int gray_value = (qRed(color) + qGreen(color) + qBlue(color)) / 3;
			ptr_image_->setPixel(i, j, qRgb(gray_value, gray_value, gray_value));
		}
	}

	update();
}

void ImageWidget::IDW() {
	p_points.clear();
	q_points.clear();
	RBF_status_ = false;
	IDW_status_ = true;
	update();
}

void ImageWidget::RBF() {
	p_points.clear();
	q_points.clear();
	IDW_status_ = false;
	RBF_status_ = true;
	update();
}

void ImageWidget::Interpolate() {
	if (p_points.empty()) {
		QMessageBox::warning(this, tr("Warning"), tr("Please select a interpolation method!"));
		return;
	}
	bool* interpolated = new bool[ptr_image_->width() * ptr_image_->height()];
	for (int i = 0; i < ptr_image_->width() * ptr_image_->height(); ++i) {
		interpolated[i] = false;
	}
	ANNpointArray dataPts = annAllocPts(ptr_image_->width() * ptr_image_->height(), 2);
	int interpolatedNum = 0;
	if (IDW_status_) {
		double mu = 2.0;
		QImage* image = new QImage(ptr_image_->width(), ptr_image_->height(), QImage::Format_RGB32);
		for (int i = 0; i < ptr_image_->width(); ++i) {
			for (int j = 0; j < ptr_image_->height(); ++j) {
				double x = 0.0, y = 0.0;
				double sigma_sum = 0.0;
				bool flag = false;
				for (int k = 0; k < p_points.size(); ++k) {
					double distance = sqrt(pow(i - p_points[k].x(), 2) + pow(j - p_points[k].y(), 2));
					if (Eigen::Vector2i(i, j) == p_points[k]) {
						flag = true;
						image->setPixel(q_points[k].x(), q_points[k].y(), ptr_image_->pixel(i, j));
						interpolated[q_points[k].x() * ptr_image_->height() + q_points[k].y()] = true;
						dataPts[interpolatedNum][0] = q_points[k].x();
						dataPts[interpolatedNum++][1] = q_points[k].y();
						break;
					}
					// D_i={0, 0, 0, 0}
					x += pow(distance, -mu) * q_points[k].x();
					y += pow(distance, -mu) * q_points[k].y();
					sigma_sum += pow(distance, -mu);
				}
				if (!flag) {
					x /= sigma_sum;
					y /= sigma_sum;
					image->setPixel((int)x, (int)y, ptr_image_->pixel(i, j));
					interpolated[(int)x * ptr_image_->height() + (int)y] = true;
					dataPts[interpolatedNum][0] = (int)x;
					dataPts[interpolatedNum++][1] = (int)y;
				}
			}
		}
		*(ptr_image_) = *image;
	}
	else if (RBF_status_) {
		double mu = 1.0;
		// 半径平方
		std::vector<double> radius_2;
		for (int i = 0; i < p_points.size(); ++i) {
			double min = 1000000;
			for (int j = 0; j < p_points.size(); ++j) {
				double distance_2 = pow(p_points[i].x() - p_points[j].x(), 2) + pow(p_points[i].y() - p_points[j].y(), 2);
				if (distance_2 < min)min = distance_2;
			}
			radius_2.push_back(min);
		}
		// q - p 矩阵
		Eigen::MatrixX2d qMinusP(p_points.size(), 2);
		for (int i = 0; i < p_points.size(); ++i) {
			qMinusP(i, 0) = double(q_points[i].x() - p_points[i].x());
			qMinusP(i, 1) = double(q_points[i].y() - p_points[i].y());
		}
		// 径向基函数
		Eigen::MatrixXd R(p_points.size(), p_points.size());
		for (int i = 0; i < p_points.size(); ++i) {
			for (int j = 0; j < p_points.size(); ++j) {
				double distance_2 = pow(p_points[i].x() - p_points[j].x(), 2) + pow(p_points[i].y() - p_points[j].y(), 2);
				R(i, j) = pow(distance_2 + radius_2[i], mu / 2);
			}
		}
		// R alpha = qMinusP
		Eigen::MatrixX2d alpha = R.colPivHouseholderQr().solve(qMinusP);
		QImage* image = new QImage(ptr_image_->width(), ptr_image_->height(), QImage::Format_RGB32);
		for (int i = 0; i < ptr_image_->width(); ++i) {
			for (int j = 0; j < ptr_image_->height(); ++j) {
				double x = i, y = j;
				bool flag = false;
				for (int k = 0; k < p_points.size(); ++k) {
					double distance_2 = pow(i - p_points[k].x(), 2) + pow(j - p_points[k].y(), 2);
					if (Eigen::Vector2i(i, j) == p_points[k]) {
						flag = true;
						image->setPixel(q_points[k].x(), q_points[k].y(), ptr_image_->pixel(i, j));
						interpolated[q_points[k].x() * ptr_image_->height() + q_points[k].y()] = true;
						dataPts[interpolatedNum][0] = q_points[k].x();
						dataPts[interpolatedNum++][1] = q_points[k].y();
						break;
					}
					x += pow(distance_2 + radius_2[k], mu / 2) * alpha(k, 0);
					y += pow(distance_2 + radius_2[k], mu / 2) * alpha(k, 1);
				}
				if (!flag) {
					if ((int)x >= 0 && (int)x < ptr_image_->width() && (int)y >= 0 && (int)y < ptr_image_->height()) {
						image->setPixel((int)x, (int)y, ptr_image_->pixel(i, j));
						interpolated[(int)x * ptr_image_->height() + (int)y] = true;
						dataPts[interpolatedNum][0] = (int)x;
						dataPts[interpolatedNum++][1] = (int)y;
					}
				}
			}
		}
		*(ptr_image_) = *image;
	}
	else {
		QMessageBox::warning(this, tr("Warning"), tr("Please select a interpolation method!"));
	}
	ANNkd_tree* kdTree = new ANNkd_tree(dataPts, interpolatedNum, 2);
	for (int i = 0; i < ptr_image_->width(); ++i) {
		for (int j = 0; j < ptr_image_->height(); ++j) {
			if (!interpolated[i * ptr_image_->height() + j]) {
				ANNpoint queryPt = annAllocPt(2);
				queryPt[0] = i;
				queryPt[1] = j;
				ANNidxArray nnIdx = new ANNidx[1];
				ANNdistArray dists = new ANNdist[1];
				kdTree->annkSearch(queryPt, 1, nnIdx, dists);
				ptr_image_->setPixel(i, j, ptr_image_->pixel(dataPts[nnIdx[0]][0], dataPts[nnIdx[0]][1]));
				delete[] nnIdx;
				delete[] dists;
				annDeallocPt(queryPt);
			}
		}
	}
	annDeallocPts(dataPts);
	delete kdTree;
	p_points.clear();
	q_points.clear();
	rect_list_.clear();
	line_list_.clear();
	IDW_status_ = false;
	RBF_status_ = false;
	update();
}

void ImageWidget::Restore() {
	*(ptr_image_) = *(ptr_image_backup_);
	p_points.clear();
	q_points.clear();
	rect_list_.clear();
	line_list_.clear();
	update();
}

void ImageWidget::mousePressEvent(QMouseEvent* event) {
	if (event->button() == Qt::LeftButton) {
		start_point_ = event->pos();
		end_point_ = event->pos();
		if (IDW_status_ || RBF_status_) {
			QPoint point = event->pos() - QPoint((width() - ptr_image_->width()) / 2, (height() - ptr_image_->height()) / 2);
			// cout << qRed(ptr_image_->pixel(point)) << ' ' << qGreen(ptr_image_->pixel(point)) << ' ' << qBlue(ptr_image_->pixel(point)) << endl;
			p_points.push_back(Eigen::Vector2i(point.x(), point.y()));
			draw_status_ = true;
		}
	}
	update();
}

void ImageWidget::mouseMoveEvent(QMouseEvent* event) {
	end_point_ = event->pos();
	update();
}

void ImageWidget::mouseReleaseEvent(QMouseEvent* event) {
	if (event->button() == Qt::LeftButton) {
		end_point_ = event->pos();
		if (IDW_status_ || RBF_status_) {
			QPoint point = event->pos() - QPoint((width() - ptr_image_->width()) / 2, (height() - ptr_image_->height()) / 2);
			q_points.push_back(Eigen::Vector2i(point.x(), point.y()));
			draw_status_ = false;
			rect_list_.push_back(QRect(start_point_.x() - 5, start_point_.y() - 5, 10, 10));
			rect_list_.push_back(QRect(end_point_.x() - 5, end_point_.y() - 5, 10, 10));
			line_list_.push_back(QLine(start_point_, end_point_));
		}
	}
	update();
}