//
// Created by Pawulon on 13/12/2023.
//

#include "Event.hpp"

Event::Event(int _index, int _type, const Vertex2D &_point) :
        index(_index), type(_type), point(_point), arc(nullptr) {}


EventPtr checkCircleEvent(bl::BLNodePtr n1, bl::BLNodePtr n2, bl::BLNodePtr n3,
                          const std::vector<Vertex2D> &points, double sweepline) {

    if (n1 == nullptr || n2 == nullptr || n3 == nullptr)
        return nullptr;


    Vertex2D p1 = points[n1->get_id()];
    Vertex2D p2 = points[n2->get_id()];
    Vertex2D p3 = points[n3->get_id()];
    Vertex2D center, bottom;

    if (p2.y > p1.y && p2.y > p3.y)
        return nullptr;

    if (!findCircleCenter(p1, p2, p3, center))
        return nullptr;

    bottom = center;
    bottom.y += (center - p2).norm();

    // check circle event
    if (fabs(bottom.y - sweepline) < POINT_EPSILON || sweepline < bottom.y) {
        // create a circle event structure
        EventPtr e = std::make_shared<Event>(-1, Event::CIRCLE, bottom);
        // initialize attributes
        e->vertexIindices = {n1->get_id(), n2->get_id(), n3->get_id()};
        e->center = center;
        e->arc = n2;
        // add reference in the corresponding node
        n2->circle_event = e;
        return e;
    }

    return nullptr;
}