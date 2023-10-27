#pragma once

#include "Shape.h"

class Poly : public Shape {
public:
    Poly();
    ~Poly();

    void Draw(QPainter& painter);

    QVector<QPoint> points;
};
