#include "Poly.h"

Poly::Poly() {
}

Poly::~Poly() {
}

void Poly::Draw(QPainter& painter) {
    painter.drawPolygon(points);
}
