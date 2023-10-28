#include "ImageWidget.h"
#include <QImage>
#include <QPainter>
#include <QtWidgets> 
#include <iostream>

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
						break;
					}
					x += pow(distance, -mu) * q_points[k].x();
					y += pow(distance, -mu) * q_points[k].y();
					sigma_sum += pow(distance, -mu);
				}
				if (!flag) {
					x /= sigma_sum;
					y /= sigma_sum;
					image->setPixel((int)x, (int)y, ptr_image_->pixel(i, j));
				}
			}
		}
		*(ptr_image_) = *image;
	}
	else if (RBF_status_) {

	}
	else {
		QMessageBox::warning(this, tr("Warning"), tr("Please select a interpolation method!"));
	}
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
		if (IDW_status_) {
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
	if (IDW_status_) {

	}

	update();
}

void ImageWidget::mouseReleaseEvent(QMouseEvent* event) {
	if (event->button() == Qt::LeftButton) {
		end_point_ = event->pos();
		if (IDW_status_) {
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