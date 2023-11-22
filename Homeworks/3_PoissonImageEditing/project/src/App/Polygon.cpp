#include "Polygon.h"

using namespace poissonedit;

Polygon::Polygon() {
	type_ = kPolygon;
	finish = false;
	// polygon.push_back(start);
}

Polygon::~Polygon() {
}

QPolygon Polygon::get_polygon() {
	return polygon;
}

void Polygon::update(int mode) {
	switch (mode) {
	case 0:
		finish = true;
		break;
	case 1:
		polygon.push_back(start);
		break;
	case 2:
		if (polygon.size() <= 1) {
			polygon.push_back(end);
			break;
		}
		polygon.pop_back();
		polygon.push_back(end);
		break;
	case 3:
		polygon.push_back(end);
		break;
	default:
		break;
	}
}


void Polygon::Draw(QPainter& painter) {
	if (finish)
		painter.drawPolygon(polygon);
	else
		painter.drawPolyline(polygon);
}
