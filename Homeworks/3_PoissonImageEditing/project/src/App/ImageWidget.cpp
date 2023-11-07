#include "ImageWidget.h"
#include <QImage>
#include <QPainter>
#include <QtWidgets> 
#include <iostream>
#include "ChildWindow.h"

using std::cout;
using std::endl;

ImageWidget::ImageWidget(ChildWindow* relatewindow) {
	// image_ = new QImage();
	// image_backup_ = new QImage();

	draw_status_ = kNone;
	is_choosing_ = false;
	is_pasting_ = false;

	point_start_ = QPoint(0, 0);
	point_end_ = QPoint(0, 0);

	source_window_ = NULL;
}

ImageWidget::~ImageWidget(void) {
}

int ImageWidget::ImageWidth() {
	return image_.cols;
}

int ImageWidget::ImageHeight() {
	return image_.rows;
}

void ImageWidget::set_draw_status_to_choose_rect() {
	draw_status_ = kChooseRect;
}

void ImageWidget::set_draw_status_to_choose_polygon() {
	draw_status_ = kChoosePolygon;
}

void ImageWidget::set_draw_status_to_paste() {
	draw_status_ = kPaste;
}

void ImageWidget::set_draw_status_to_paste_possion() {
	draw_status_ = kPastePossion;
}

void ImageWidget::set_draw_status_to_paste_mixed() {
	draw_status_ = kPasteMixed;
}

cv::Mat ImageWidget::image() {
	return image_;
}

void ImageWidget::set_source_window(ChildWindow* childwindow) {
	source_window_ = childwindow;
}

void ImageWidget::paintEvent(QPaintEvent* paintevent) {
	QPainter painter;
	painter.begin(this);

	// Draw background
	painter.setBrush(Qt::lightGray);
	QRect back_rect(0, 0, width(), height());
	painter.drawRect(back_rect);

	// Draw image
	QRect rect = QRect(0, 0, image_.cols, image_.rows);
	QImage Qimage = QImage((unsigned char*)(image_.data), image_.cols, image_.rows, image_.step, QImage::Format_RGB888);
	painter.drawImage(rect, Qimage);

	// Draw choose region
	painter.setBrush(Qt::NoBrush);
	painter.setPen(Qt::red);
	switch (draw_status_) {
	case kChooseRect:
		painter.drawRect(point_start_.x(), point_start_.y(),
			point_end_.x() - point_start_.x(), point_end_.y() - point_start_.y());
		break;
	case kChoosePolygon:
		painter.drawPolygon(polygon_points_);
		break;
	}

	painter.end();
}

void ImageWidget::mousePressEvent(QMouseEvent* mouseevent) {
	if (Qt::LeftButton == mouseevent->button()) {
		switch (draw_status_) {
		case kChooseRect:
			is_choosing_ = true;
			point_start_ = point_end_ = mouseevent->pos();
			break;

		case kChoosePolygon:
			is_choosing_ = true;
			polygon_points_.clear();
			polygon_points_.push_back(mouseevent->pos());
			break;

		case kPaste:
		{
			is_pasting_ = true;

			switch (source_window_->imagewidget_->draw_status_) {
			case kChooseRect:
			{
				// Start point in object image
				int xpos = mouseevent->pos().rx();
				int ypos = mouseevent->pos().ry();

				// Start point in source image
				int xsourcepos = source_window_->imagewidget_->point_start_.rx();
				int ysourcepos = source_window_->imagewidget_->point_start_.ry();

				// Width and Height of rectangle region
				int w = source_window_->imagewidget_->point_end_.rx()
					- source_window_->imagewidget_->point_start_.rx() + 1;
				int h = source_window_->imagewidget_->point_end_.ry()
					- source_window_->imagewidget_->point_start_.ry() + 1;

				// Paste

				// Paste
				for (int i = 0; i < w; i++) {
					for (int j = 0; j < h; j++) {
						if (xpos + i >= 0 && xpos + i < image_.cols && ypos + j >= 0 && ypos + j < image_.rows)
							image_.at<cv::Vec3b>(ypos + j, xpos + i) = source_window_->imagewidget_->image().at<cv::Vec3b>(ysourcepos + j, xsourcepos + i);
						// image_->setPixel(xpos + i, ypos + j, source_window_->imagewidget_->image()->pixel(xsourcepos + i, ysourcepos + j));
					}
				}

				break;
			}
			case kChoosePolygon:
			{
				// Start point in object image
				int xpos = mouseevent->pos().rx();
				int ypos = mouseevent->pos().ry();

				// Points of polygon region 扫描线算法
				// 求多边形包围盒
				int xsourcemin = source_window_->imagewidget_->polygon_points_[0].rx();
				int ysourcemin = source_window_->imagewidget_->polygon_points_[0].ry();
				int xsourcemax = source_window_->imagewidget_->polygon_points_[0].rx();
				int ysourcemax = source_window_->imagewidget_->polygon_points_[0].ry();
				for (int i = 0; i < source_window_->imagewidget_->polygon_points_.size(); ++i) {
					if (source_window_->imagewidget_->polygon_points_[i].rx() < xsourcemin) {
						xsourcemin = source_window_->imagewidget_->polygon_points_[i].rx();
					}
					if (source_window_->imagewidget_->polygon_points_[i].ry() < ysourcemin) {
						ysourcemin = source_window_->imagewidget_->polygon_points_[i].ry();
					}
					if (source_window_->imagewidget_->polygon_points_[i].rx() > xsourcemax) {
						xsourcemax = source_window_->imagewidget_->polygon_points_[i].rx();
					}
					if (source_window_->imagewidget_->polygon_points_[i].ry() > ysourcemax) {
						ysourcemax = source_window_->imagewidget_->polygon_points_[i].ry();
					}
				}
				// 扫描线
				// 边表
				struct Edge {
					double x;			// 下端点x坐标
					int ymax;			// 上端点y坐标
					double dx;			// dx/dy
				};
				std::vector<std::list<Edge>> edge_table(ysourcemax - ysourcemin + 1);	// 按下端点的y坐标排序，相同的按x坐标排序
				// 建立边表
				for (int i = 0; i < source_window_->imagewidget_->polygon_points_.size(); ++i) {
					Edge* edge = new Edge;
					QPoint p1 = source_window_->imagewidget_->polygon_points_[i];
					QPoint p2 = source_window_->imagewidget_->polygon_points_[(i + 1) % source_window_->imagewidget_->polygon_points_.size()];
					if (p1.ry() < p2.ry()) {
						edge->x = p1.rx();
						edge->ymax = p2.ry();
						edge->dx = (double)(p2.rx() - p1.rx()) / (double)(p2.ry() - p1.ry());
						auto tmp = edge_table[p1.ry() - ysourcemin].begin();
						while (tmp != edge_table[p1.ry() - ysourcemin].end() && tmp->x < edge->x) {
							++tmp;
						}
						edge_table[p1.ry() - ysourcemin].insert(tmp, *edge);
					}
					else {
						edge->x = p2.rx();
						edge->ymax = p1.ry();
						edge->dx = (double)(p1.rx() - p2.rx()) / (double)(p1.ry() - p2.ry());
						auto tmp = edge_table[p2.ry() - ysourcemin].begin();
						while (tmp != edge_table[p2.ry() - ysourcemin].end() && tmp->x < edge->x) {
							++tmp;
						}
						edge_table[p2.ry() - ysourcemin].insert(tmp, *edge);
					}
				}
				// 扫描线填充
				// 建立活动边表
				std::list<Edge> active_edge_table;
				auto tmp = edge_table[0].begin();
				while (tmp != edge_table[0].end()) {
					active_edge_table.push_back(*tmp);
					++tmp;
				}
				for (int i = 0; i < edge_table.size(); ++i) {
					// 填充
					tmp = active_edge_table.begin();
					while (tmp != active_edge_table.end()) {
						int x1 = tmp->x;
						int x2 = (++tmp)->x;
						for (int j = x1; j <= x2; ++j) {
							if (xpos + j - xsourcemin >= 0 && xpos + j - xsourcemin < image_.cols && ypos + i >= 0 && ypos + i < image_.rows)
								image_.at<cv::Vec3b>(ypos + i, xpos + j - xsourcemin) = source_window_->imagewidget_->image().at<cv::Vec3b>(i + ysourcemin, j);
							// image_->setPixel(xpos + j - xsourcemin, ypos + i, source_window_->imagewidget_->image()->pixel(j, i + ysourcemin));
						}
						++tmp;
					}
					// 删去活动边表中不再相交的边，更新活动边表中的x
					tmp = active_edge_table.begin();
					while (tmp != active_edge_table.end()) {
						if (tmp->ymax == i + 1 + ysourcemin)
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
				break;
			}
			}
			break;
		}

		case kPastePossion:
		{
			//TODO
			is_pasting_ = true;
			switch (source_window_->imagewidget_->draw_status_) {
			case kChooseRect:
			{
				int xmin = std::min(source_window_->imagewidget_->point_start_.rx(), source_window_->imagewidget_->point_end_.rx());
				int ymin = std::min(source_window_->imagewidget_->point_start_.ry(), source_window_->imagewidget_->point_end_.ry());
				int xmax = std::max(source_window_->imagewidget_->point_start_.rx(), source_window_->imagewidget_->point_end_.rx());
				int ymax = std::max(source_window_->imagewidget_->point_start_.ry(), source_window_->imagewidget_->point_end_.ry());
				Eigen::Vector2i rect_start(xmin, ymin);
				Eigen::Vector2i rect_end(xmax, ymax);
				PoissonSolver solver(rect_start, rect_end, source_window_->imagewidget_->image(), image_, mouseevent->pos().rx(), mouseevent->pos().ry());
				image_ = solver.solve();
				break;
			}
			case kChoosePolygon:
			{

			}

			}
			break;
		}

		case kPasteMixed:
		{
			//TODO
			is_pasting_ = true;
			switch (source_window_->imagewidget_->draw_status_) {
			case kChooseRect:
			{

			}
			case kChoosePolygon:
			{

			}

			}
			break;
		}

		default:
			break;
		}
	}

	if (Qt::RightButton == mouseevent->button()) {
		if (draw_status_ == kChoosePolygon) {
			if (is_choosing_) {
				polygon_points_.push_back(mouseevent->pos());
			}
		}
	}
	update();
}

void ImageWidget::mouseMoveEvent(QMouseEvent* mouseevent) {
	switch (draw_status_) {
	case kChooseRect:
		// Store point position for rectangle region
		if (is_choosing_) {
			point_end_ = mouseevent->pos();
		}
		break;

	case kChoosePolygon:
		// Store point position for polygon region
		if (is_choosing_) {
			if (polygon_points_.size() == 1)
				polygon_points_.push_back(mouseevent->pos());
			polygon_points_.pop_back();
			polygon_points_.push_back(mouseevent->pos());
		}
		break;

	case kPaste:
		// Paste rectangle region to object image
		if (is_pasting_) {
			if (source_window_->imagewidget_->draw_status_ == kChooseRect) {
				// Start point in object image
				int xpos = mouseevent->pos().rx();
				int ypos = mouseevent->pos().ry();

				// Start point in source image
				int xsourcepos = source_window_->imagewidget_->point_start_.rx();
				int ysourcepos = source_window_->imagewidget_->point_start_.ry();

				// Width and Height of rectangle region
				int w = source_window_->imagewidget_->point_end_.rx()
					- source_window_->imagewidget_->point_start_.rx() + 1;
				int h = source_window_->imagewidget_->point_end_.ry()
					- source_window_->imagewidget_->point_start_.ry() + 1;

				// Paste

				// Restore image 
				image_ = image_backup_.clone();

				// Paste
				for (int i = 0; i < w; i++) {
					for (int j = 0; j < h; j++) {
						if (xpos + i >= 0 && xpos + i < image_.cols && ypos + j >= 0 && ypos + j < image_.rows)
							image_.at<cv::Vec3b>(ypos + j, xpos + i) = source_window_->imagewidget_->image().at<cv::Vec3b>(ysourcepos + j, xsourcepos + i);
						// image_->setPixel(xpos + i, ypos + j, source_window_->imagewidget_->image()->pixel(xsourcepos + i, ysourcepos + j));
					}
				}

			}
			if (source_window_->imagewidget_->draw_status_ == kChoosePolygon) {
				// Start point in object image
				int xpos = mouseevent->pos().rx();
				int ypos = mouseevent->pos().ry();

				// Points of polygon region 扫描线算法
				// 求多边形包围盒
				int xsourcemin = source_window_->imagewidget_->polygon_points_[0].rx();
				int ysourcemin = source_window_->imagewidget_->polygon_points_[0].ry();
				int xsourcemax = source_window_->imagewidget_->polygon_points_[0].rx();
				int ysourcemax = source_window_->imagewidget_->polygon_points_[0].ry();
				for (int i = 0; i < source_window_->imagewidget_->polygon_points_.size(); ++i) {
					if (source_window_->imagewidget_->polygon_points_[i].rx() < xsourcemin) {
						xsourcemin = source_window_->imagewidget_->polygon_points_[i].rx();
					}
					if (source_window_->imagewidget_->polygon_points_[i].ry() < ysourcemin) {
						ysourcemin = source_window_->imagewidget_->polygon_points_[i].ry();
					}
					if (source_window_->imagewidget_->polygon_points_[i].rx() > xsourcemax) {
						xsourcemax = source_window_->imagewidget_->polygon_points_[i].rx();
					}
					if (source_window_->imagewidget_->polygon_points_[i].ry() > ysourcemax) {
						ysourcemax = source_window_->imagewidget_->polygon_points_[i].ry();
					}
				}
				// 扫描线
				// 边表
				struct Edge {
					double x;			// 下端点x坐标
					int ymax;			// 上端点y坐标
					double dx;			// dx/dy
				};
				std::vector<std::list<Edge>> edge_table(ysourcemax - ysourcemin + 1);	// 按下端点的y坐标排序，相同的按x坐标排序
				// 建立边表
				for (int i = 0; i < source_window_->imagewidget_->polygon_points_.size(); ++i) {
					Edge* edge = new Edge;
					QPoint p1 = source_window_->imagewidget_->polygon_points_[i];
					QPoint p2 = source_window_->imagewidget_->polygon_points_[(i + 1) % source_window_->imagewidget_->polygon_points_.size()];
					if (p1.ry() < p2.ry()) {
						edge->x = p1.rx();
						edge->ymax = p2.ry();
						edge->dx = (double)(p2.rx() - p1.rx()) / (double)(p2.ry() - p1.ry());
						auto tmp = edge_table[p1.ry() - ysourcemin].begin();
						while (tmp != edge_table[p1.ry() - ysourcemin].end() && tmp->x < edge->x) {
							++tmp;
						}
						edge_table[p1.ry() - ysourcemin].insert(tmp, *edge);
					}
					else {
						edge->x = p2.rx();
						edge->ymax = p1.ry();
						edge->dx = (double)(p1.rx() - p2.rx()) / (double)(p1.ry() - p2.ry());
						auto tmp = edge_table[p2.ry() - ysourcemin].begin();
						while (tmp != edge_table[p2.ry() - ysourcemin].end() && tmp->x < edge->x) {
							++tmp;
						}
						edge_table[p2.ry() - ysourcemin].insert(tmp, *edge);
					}
				}
				// 扫描线填充
				// 建立活动边表
				std::list<Edge> active_edge_table;
				auto tmp = edge_table[0].begin();
				while (tmp != edge_table[0].end()) {
					active_edge_table.push_back(*tmp);
					++tmp;
				}
				image_ = image_backup_.clone();
				for (int i = 0; i < edge_table.size(); ++i) {
					// 填充
					tmp = active_edge_table.begin();
					while (tmp != active_edge_table.end()) {
						int x1 = tmp->x;
						int x2 = (++tmp)->x;
						for (int j = x1; j <= x2; ++j) {
							if (xpos + j - xsourcemin >= 0 && xpos + j - xsourcemin < image_.cols && ypos + i >= 0 && ypos + i < image_.rows)
								image_.at<cv::Vec3b>(ypos + i, xpos + j - xsourcemin) = source_window_->imagewidget_->image().at<cv::Vec3b>(i + ysourcemin, j);
							// image_->setPixel(xpos + j - xsourcemin, ypos + i, source_window_->imagewidget_->image()->pixel(j, i + ysourcemin));
						}
						++tmp;
					}
					// 删去活动边表中不再相交的边，更新活动边表中的x
					tmp = active_edge_table.begin();
					while (tmp != active_edge_table.end()) {
						if (tmp->ymax == i + 1 + ysourcemin)
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
		}
		break;

	default:
		break;
	}

	update();
}

void ImageWidget::mouseReleaseEvent(QMouseEvent* mouseevent) {
	if (Qt::LeftButton == mouseevent->button()) {
		switch (draw_status_) {
		case kChooseRect:
			if (is_choosing_) {
				point_end_ = mouseevent->pos();
				is_choosing_ = false;
			}
			break;

		case kChoosePolygon:
			if (is_choosing_) {
				is_choosing_ = false;
				polygon_points_.pop_back();
			}
			break;

		case kPaste:
			if (is_pasting_) {
				is_pasting_ = false;
			}
			break;

		case kPasteMixed:
			if (is_pasting_) {
				is_pasting_ = false;
			}
			break;

		default:
			break;
		}
	}
	update();
}

void ImageWidget::Open(QString filename) {
	// Load file
	if (!filename.isEmpty()) {
		image_ = cv::imread(filename.toLatin1().data());
		cv::cvtColor(image_, image_, cv::COLOR_BGR2RGB);
		image_backup_ = image_.clone();
	}

	//	setFixedSize(image_->width(), image_->height());
	//	relate_window_->setWindowFlags(Qt::Dialog);
	//	relate_window_->setFixedSize(QSize(image_->width(), image_->height()));
	//	relate_window_->setWindowFlags(Qt::SubWindow);

		//image_->invertPixels(QImage::InvertRgb);
		//*(image_) = image_->mirrored(true, true);
		//*(image_) = image_->rgbSwapped();
	cout << "image size: " << image_.cols << ' ' << image_.rows << endl;
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

	cv::imwrite(filename.toLatin1().data(), image_);
}

void ImageWidget::Invert() {
	for (int i = 0; i < image_.cols; i++) {
		for (int j = 0; j < image_.rows; j++) {
			// QRgb color = image_->pixel(i, j);
			// image_->setPixel(i, j, qRgb(255 - qRed(color), 255 - qGreen(color), 255 - qBlue(color)));
			image_.at<cv::Vec3b>(j, i)[0] = 255 - image_.at<cv::Vec3b>(j, i)[0];
			image_.at<cv::Vec3b>(j, i)[1] = 255 - image_.at<cv::Vec3b>(j, i)[1];
			image_.at<cv::Vec3b>(j, i)[2] = 255 - image_.at<cv::Vec3b>(j, i)[2];
		}
	}

	// equivalent member function of class QImage
	// image_->invertPixels(QImage::InvertRgb);
	update();
}

void ImageWidget::Mirror(bool ishorizontal, bool isvertical) {
	// QImage image_tmp(*(image_));
	cv::Mat image_tmp = image_.clone();
	int cols = image_.cols;
	int rows = image_.rows;

	if (ishorizontal) {
		if (isvertical) {
			for (int i = 0; i < rows; i++) {
				for (int j = 0; j < cols; j++) {
					image_.at<cv::Vec3b>(i, j)[0] = image_tmp.at<cv::Vec3b>(rows - 1 - i, cols - 1 - j)[0];
					image_.at<cv::Vec3b>(i, j)[1] = image_tmp.at<cv::Vec3b>(rows - 1 - i, cols - 1 - j)[1];
					image_.at<cv::Vec3b>(i, j)[2] = image_tmp.at<cv::Vec3b>(rows - 1 - i, cols - 1 - j)[2];
					// image_->setPixel(i, j, image_tmp.pixel(width - 1 - i, height - 1 - j));
				}
			}
		}
		else {
			for (int i = 0; i < rows; i++) {
				for (int j = 0; j < cols; j++) {
					image_.at<cv::Vec3b>(i, j)[0] = image_tmp.at<cv::Vec3b>(rows - 1 - i, j)[0];
					image_.at<cv::Vec3b>(i, j)[1] = image_tmp.at<cv::Vec3b>(rows - 1 - i, j)[1];
					image_.at<cv::Vec3b>(i, j)[2] = image_tmp.at<cv::Vec3b>(rows - 1 - i, j)[2];
					// image_->setPixel(i, j, image_tmp.pixel(i, height - 1 - j));
				}
			}
		}

	}
	else {
		if (isvertical) {
			for (int i = 0; i < rows; i++) {
				for (int j = 0; j < cols; j++) {
					image_.at<cv::Vec3b>(i, j)[0] = image_tmp.at<cv::Vec3b>(i, cols - 1 - j)[0];
					image_.at<cv::Vec3b>(i, j)[1] = image_tmp.at<cv::Vec3b>(i, cols - 1 - j)[1];
					image_.at<cv::Vec3b>(i, j)[2] = image_tmp.at<cv::Vec3b>(i, cols - 1 - j)[2];
					// image_->setPixel(i, j, image_tmp.pixel(width - 1 - i, j));
				}
			}
		}
	}

	// equivalent member function of class QImage
	//*(image_) = image_->mirrored(true, true);
	update();
}

void ImageWidget::TurnGray() {
	for (int i = 0; i < image_.cols; i++) {
		for (int j = 0; j < image_.rows; j++) {
			// QRgb color = image_->pixel(i, j);
			// int gray_value = (qRed(color) + qGreen(color) + qBlue(color)) / 3;
			// image_->setPixel(i, j, qRgb(gray_value, gray_value, gray_value));
			int gray_value = (image_.at<cv::Vec3b>(j, i)[0] + image_.at<cv::Vec3b>(j, i)[1] + image_.at<cv::Vec3b>(j, i)[2]) / 3;
			image_.at<cv::Vec3b>(j, i)[0] = gray_value;
			image_.at<cv::Vec3b>(j, i)[1] = gray_value;
			image_.at<cv::Vec3b>(j, i)[2] = gray_value;
		}
	}

	update();
}

void ImageWidget::Restore() {
	image_ = image_backup_.clone();
	point_start_ = point_end_ = QPoint(0, 0);
	update();
}
