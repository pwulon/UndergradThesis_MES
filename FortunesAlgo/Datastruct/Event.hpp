//
// Created by Pawulon on 13/12/2023.
//

#ifndef MES_EVENT_HPP
#define MES_EVENT_HPP

#include "Beachline.hpp"
#include "../Math/Circle.hpp"

namespace mes::fortunes{
namespace bl = beachline;


struct Event {
    enum {
        SITE = 0, CIRCLE = 1, SKIP = 2,
    };

    int type;
    Vertex2D point;

    /*
     Site event attributes:
     */
    int index;
    std::vector<int> vertexIindices;

    /*
     Circle event attributes:
     */
    Vertex2D center;
    bl::BLNodePtr arc;


    Event(int _index = -1, int _type = Event::SKIP, const Vertex2D &_point = Vertex2D(0.0, 0.0));

};


typedef std::shared_ptr<Event> EventPtr;


struct Point2DComparator {
    bool operator()(const Vertex2D &p1, const Vertex2D &p2) {
        return (p1.y == p2.y && p1.x > p2.x) || p1.y > p2.y;
    }
};

struct Point2DComparator2 {
    bool operator()(const Vertex2D &p1, const Vertex2D &p2) {
        return (p1.y == p2.y && p1.x < p2.x) || p1.y < p2.y;
    }
};


struct EventPtrComparator {
    Point2DComparator point_cmp;

    bool operator()(const EventPtr &e1, const EventPtr &e2) {
        return point_cmp(e1->point, e2->point);
    }
};


EventPtr checkCircleEvent(bl::BLNodePtr n1, bl::BLNodePtr n2, bl::BLNodePtr n3,
                          const std::vector<Vertex2D> &points, double sweepline);

}
#endif //MES_EVENT_HPP
