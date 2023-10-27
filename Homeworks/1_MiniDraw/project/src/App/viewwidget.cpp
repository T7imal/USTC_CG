#include "viewwidget.h"

ViewWidget::ViewWidget(QWidget* parent)
	: QWidget(parent) {
	ui.setupUi(this);
	draw_status_ = false;
	shape_ = NULL;
	type_ = Shape::kDefault;
}

ViewWidget::~ViewWidget() {
}

void ViewWidget::setLine() {
	type_ = Shape::kLine;
}

void ViewWidget::setRect() {
	type_ = Shape::kRect;
}

void ViewWidget::setElli() {
	type_ = Shape::kElli;
}

void ViewWidget::setPoly() {
	type_ = Shape::kPoly;
}

void ViewWidget::mousePressEvent(QMouseEvent* event) {
	if (Qt::LeftButton == event->button()) {
		switch (type_) {
		case Shape::kLine:
			shape_ = new Line();
			break;
		case Shape::kDefault:
			// shape_保持原狀
			break;

		case Shape::kRect:
			shape_ = new Rect();
			break;

		case Shape::kElli:
			shape_ = new Elli();
			break;

		case Shape::kPoly:
			shape_ = new Poly();
			Poly* poly = dynamic_cast<Poly*>(shape_);
			poly->points.push_back(event->pos());
			break;
		}
		if (shape_ != NULL) {
			draw_status_ = true;
			start_point_ = end_point_ = event->pos();
			shape_->set_start(start_point_);
			shape_->set_end(end_point_);
		}
	}
	if (Qt::RightButton == event->button()) {
		if (type_ == Shape::kPoly && shape_ != NULL) {
			Poly* poly = dynamic_cast<Poly*>(shape_);
			poly->points.push_back(event->pos());
		}
	}
	update();
}

void ViewWidget::mouseMoveEvent(QMouseEvent* event) {
	if (draw_status_ && shape_ != NULL) {
		end_point_ = event->pos();
		shape_->set_end(end_point_);
	}
	if (type_ == Shape::kPoly && shape_ != NULL) {
		Poly* poly = dynamic_cast<Poly*>(shape_);
		if (poly->points.size() == 1)
			poly->points.push_back(event->pos());
		poly->points.pop_back();
		poly->points.push_back(event->pos());
	}
}

void ViewWidget::mouseReleaseEvent(QMouseEvent* event) {
	if (Qt::LeftButton == event->button() && shape_ != NULL) {
		draw_status_ = false;
		shape_list_.push_back(shape_);
		shape_ = NULL;
	}
}

void ViewWidget::paintEvent(QPaintEvent*) {
	QPainter painter(this);

	for (int i = 0; i < shape_list_.size(); i++) {
		shape_list_[i]->Draw(painter);
	}

	if (shape_ != NULL) {
		shape_->Draw(painter);
	}

	update();
}