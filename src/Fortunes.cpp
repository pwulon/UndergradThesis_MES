//
// Created by Pawulon on 13/12/2023.
//
#include "../FortunesAlgo/Fortunes.hpp"


#define BREAKPOINTS_EPSILON 1.0e-5




void mes::fortunes::pointsRot(std::vector<Vertex2D> &points, double ang){
    double alfa = ang*M_PI/180.;
    double sin1 = sin(alfa);
    double cos1 = cos(alfa);

    for (int i = 0; i < points.size(); i++) {
        double x = points[i].x*cos1 - points[i].y*sin1;
        double y = points[i].x*sin1 + points[i].y*cos1;
        points[i].x = x;
        points[i].y = y;
    }
}

void mes::fortunes::build(std::vector<Vertex2D> &points,
           std::vector<mes::ElementIndices> &elements,
           std::vector<Wall> walls, mes::baseFuncType _bft, bool withPointRot) {

    // rotate points to doge edge case when too many points are in a straight line
    if(withPointRot) {
        pointsRot(points, 1);
        for(auto &w:walls) w.rot(1);
    }

    // create a priority queue
    std::priority_queue<EventPtr, std::vector<EventPtr>, EventPtrComparator> pq;

    // initialize it with all site events
    for (size_t i = 0; i < points.size(); ++i) {
        pq.push(std::make_shared<Event>(static_cast<int>(i), Event::SITE, points[i]));
    }

    // create a beachline tree
    bl::BLNodePtr root;
    double sweepline = 0L; // current position of the sweepline

    // process events
    while (!pq.empty()) {

        // extract new event from the queue
        EventPtr e = pq.top(); pq.pop();

        // set position of a sweepline
        sweepline = e->point.y;

        if (e->type == Event::SITE) { // handle site event

            int point_i = e->index;
            if (root == nullptr) { // init empty beachline tree
                root = std::make_shared<bl::BLNode>(std::make_pair(point_i, point_i), &sweepline, &points);

            } else { // if it's not empty
                bl::BLNodePtr arc = bl::find(root, e->point.x);
                bl::BLNodePtr subtree, left_leaf, right_leaf;

                if (arc->circle_event != nullptr) {
                    EventPtr circle_e = arc->circle_event;
                    circle_e->type = Event::SKIP; // ignore corresponding event
                }

                // check number of intersection points
                int isp_num = intersectionPointsNum(points[arc->get_id()], e->point, sweepline);

                // different subtrees depending on the number of intersection points
                if (isp_num == 1) {
                    subtree = bl::make_simple_subtree(point_i, arc->get_id(), &sweepline, &points);
                    left_leaf = subtree->left;
                    right_leaf = subtree->right;
                } else if (isp_num == 2) {
                    subtree = bl::make_subtree(point_i, arc->get_id(), &sweepline, &points);
                    left_leaf = subtree->left;
                    right_leaf = subtree->right->right;
                } else {
                    continue;
                }

                if (arc->prev != nullptr)
                    bl::connect(arc->prev, left_leaf);

                if (arc->next != nullptr)
                    bl::connect(right_leaf, arc->next);

                // Replace old leaf with a subtree and rebalance it
                root = bl::replace(arc, subtree);

                // Check circle events
                EventPtr circle_event = checkCircleEvent(left_leaf->prev, left_leaf, left_leaf->next, points, sweepline);
                if (circle_event != nullptr) {
                    pq.push(circle_event);
                }
                circle_event = checkCircleEvent(right_leaf->prev, right_leaf, right_leaf->next, points, sweepline);
                if (circle_event != nullptr) {
                    pq.push(circle_event);
                }
            }

        } else if (e->type == Event::CIRCLE) { // handle circle event

            bl::BLNodePtr arc = e->arc, prev_leaf, next_leaf;

            // get breakpoint nodes
            std::pair<bl::BLNodePtr, bl::BLNodePtr> breakpoints = bl::breakpoints(arc);

            // recheck if it's a false alarm 1
            if (breakpoints.first == nullptr || breakpoints.second == nullptr) {
                continue;
            }

            // recheck if it's a false alarm 2
            double v1 = breakpoints.first->value(), v2 = breakpoints.second->value();

            if (fabs(v1 - v2) > BREAKPOINTS_EPSILON) {
                continue;
            }

            bool inWall = false;
            for(auto& w: walls){
                if(w.isInsideWall(e->center)){
                    inWall = true;
                    elements.emplace_back(e->vertexIindices, w.type, _bft);
                    break;
                }
            }
            if(!inWall){
                elements.emplace_back(e->vertexIindices, mes::AIR, _bft);
            }



            // remove circle event corresponding to next leaf
            if (arc->prev != nullptr && arc->prev->circle_event != nullptr) {
                EventPtr circle_e = arc->prev->circle_event;
                circle_e->type = Event::SKIP; // ignore corresponding event
            }

            // remove circle event corresponding to prev leaf
            if (arc->next != nullptr && arc->next->circle_event != nullptr) {
                EventPtr circle_e = arc->next->circle_event;
                circle_e->type = Event::SKIP; // ignore corresponding event
            }

            // store pointers to the next and previous leaves
            prev_leaf = arc->prev;
            next_leaf = arc->next;

            // They should not be null
            assert(prev_leaf != nullptr);
            assert(next_leaf != nullptr);

            // get node associated with a new edge
            bl::BLNodePtr new_edge_node;
            if (arc->parent == breakpoints.first)
                new_edge_node = breakpoints.second;
            else
                new_edge_node = breakpoints.first;

            // remove arc from the beachline
            root = bl::remove(arc);

            // check new circle events
            if (prev_leaf != nullptr && next_leaf != nullptr) {
                EventPtr circle_event = checkCircleEvent(prev_leaf->prev, prev_leaf, next_leaf, points, sweepline);
                if (circle_event != nullptr) {
                    pq.push(circle_event);
                }
                circle_event = checkCircleEvent(prev_leaf, next_leaf, next_leaf->next, points, sweepline);
                if (circle_event != nullptr) {
                    pq.push(circle_event);
                }
            }
        }
    }
    if(withPointRot)  pointsRot(points, -1);
}
