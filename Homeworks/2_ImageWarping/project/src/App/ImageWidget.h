#pragma once
#include <QWidget>
#include "Eigen/Dense"

QT_BEGIN_NAMESPACE
class QImage;
class QPainter;
QT_END_NAMESPACE

class ImageWidget :
	public QWidget {
	Q_OBJECT

public:
	ImageWidget(void);
	~ImageWidget(void);

protected:
	void paintEvent(QPaintEvent* paintevent);

public slots:
	// File IO
	void Open();												// Open an image file, support ".bmp, .png, .jpg" format
	void Save();												// Save image to current file
	void SaveAs();												// Save image to another file

	// Image processing
	void Invert();												// Invert pixel value in image
	void Mirror(bool horizontal = false, bool vertical = true);		// Mirror image vertically or horizontally
	void TurnGray();											// Turn image to gray-scale map
	void IDW();													// Inverse distance-weighted interpolation methods
	void RBF();													// Radial basis functions interpolation method
	void Interpolate();
	void Restore();												// Restore image to origin

public:
	void mousePressEvent(QMouseEvent* event);
	void mouseMoveEvent(QMouseEvent* event);
	void mouseReleaseEvent(QMouseEvent* event);

private:
	QImage* ptr_image_;				// image 
	QImage* ptr_image_backup_;
	bool draw_status_ = false;
	QPoint start_point_;
	QPoint end_point_;
	std::vector<QRect> rect_list_;
	std::vector<QLine> line_list_;
	bool IDW_status_ = false;
	bool RBF_status_ = false;
	std::vector<Eigen::Vector2i> p_points;
	std::vector<Eigen::Vector2i> q_points;
};

