/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
 * Department of Computer Science.                                                 *
 * For a list of authors please refer to the file "CREDITS.txt".                   *
 *                                                                                 *
 * This file is part of the Voreen software package. Voreen is free software:      *
 * you can redistribute it and/or modify it under the terms of the GNU General     *
 * Public License version 2 as published by the Free Software Foundation.          *
 *                                                                                 *
 * Voreen is distributed in the hope that it will be useful, but WITHOUT ANY       *
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR   *
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.      *
 *                                                                                 *
 * You should have received a copy of the GNU General Public License in the file   *
 * "LICENSE.txt" along with this file. If not, see <http://www.gnu.org/licenses/>. *
 *                                                                                 *
 * For non-commercial academic use see the license exception specified in the file *
 * "LICENSE-academic.txt". To get information about commercial licensing please    *
 * contact the authors.                                                            *
 *                                                                                 *
 ***********************************************************************************/

#include "geometryclippingwidget.h"
#include "voreen/core/utils/voreenqualitymode.h"

namespace voreen {

const float GeometryClippingWidget::EPS(0.000001f);

GeometryClippingWidget::GeometryClippingWidget()
    : GeometryRendererBase()
    , inport_(Port::INPORT, "geometry", "Geometry Input")
    , renderGeometry_("renderGeometry", "Enabled", true)
    , planeNormal_("planeNormal", "Plane Normal", tgt::vec3(0, 1, 0), tgt::vec3(-1), tgt::vec3(1))
    , planePosition_("planePosition", "Plane Position", 0.0f, -10000.0f, 10000.0f, Processor::INVALID_RESULT, NumericProperty<float>::STATIC)
    , lightPosition_("lightPosition", "Light source position", tgt::vec4(2.3f, 1.5f, 1.5f, 1.f),
                     tgt::vec4(-10), tgt::vec4(10))
    , planeColor_("planeColor", "Plane Color", tgt::vec4(0.0f, 0.0f, 1.0f, 0.05f))
    , anchorColor_("anchorColor", "Anchor Color", tgt::vec4(0.0f, 0.0f, 1.0f, 1.0f))
    , grabbedAnchorColor_("grabbedAnchorColor", "Grabbed Anchor Color", tgt::vec4(1.0f, 0.0f, 0.0f, 1.0f))
    , lineColor_("lineColor", "Line Color", tgt::vec4(0.0f, 0.5f, 1.0f, 1.0f))
    , width_("lineWidth", "Line Width", 2.0f, 1.0f, 10.0f)
    , actualPlane_(-tgt::vec4(planeNormal_.get(), planePosition_.get()))
    , isAnchorGrabbed_(false)
    , crossNewAnchorFirst_(true)
    , isLineGrabbed_(false)
{
    addPort(inport_);

    moveEventProp_ = new EventProperty<GeometryClippingWidget>(
        "mouseEvent.clipplaneManipulation", "Clipplane manipulation", this,
        &GeometryClippingWidget::planeManipulation,
        tgt::MouseEvent::MOUSE_BUTTON_LEFT,
        tgt::MouseEvent::MOTION | tgt::MouseEvent::PRESSED | tgt::MouseEvent::RELEASED,
        tgt::Event::MODIFIER_NONE, false);
    addEventProperty(moveEventProp_);

    addProperty(renderGeometry_);
    addProperty(planeNormal_);
    addProperty(planePosition_);
    addProperty(planeColor_);
    addProperty(anchorColor_);
    addProperty(grabbedAnchorColor_);
    addProperty(lineColor_);
    addProperty(width_);
    addProperty(lightPosition_);
}

GeometryClippingWidget::~GeometryClippingWidget() {
    delete moveEventProp_;
}

void GeometryClippingWidget::initialize() {
    GeometryRendererBase::initialize();
}

Processor* GeometryClippingWidget::create() const {
    return new GeometryClippingWidget();
}

void GeometryClippingWidget::setIDManager(IDManager* idm) {
    if (idManager_ == idm)
        return;

    idManager_ = idm;
    if (idManager_) {
        for (std::vector<Anchor>::iterator it = anchors_.begin(); it != anchors_.end(); ++it)
            idManager_->registerObject(&(*it));
        for (std::vector<Line>::iterator it = lines_.begin(); it != lines_.end(); ++it)
            idManager_->registerObject(&(*it));
    }
}

void GeometryClippingWidget::planeManipulation(tgt::MouseEvent* e) {
    tgtAssert(inport_.getData(), "no input geometry");

    //no plane manipulation if rendering is not enabled
    if (!renderGeometry_.get())
        return;

    e->ignore();
    // Mouse button pressend (indicates clipping plane manipulation start)?
    if (e->action() & tgt::MouseEvent::PRESSED) {
        const void* obj = idManager_->getObjectAtPos(tgt::ivec2(e->coord().x, e->viewport().y - e->coord().y));

        // Check if an anchor was selected...
        for (std::vector<Anchor>::iterator it = anchors_.begin(); it != anchors_.end(); ++it)
            if (&(*it) == obj) {
                isAnchorGrabbed_ = true;
                grabbedAnchor_ = *it;

                // Store also the first two anchors which are not equal to the grabbed anchor...
                bool isSecondAnchorFound = false;
                for (std::vector<Anchor>::iterator anchor = anchors_.begin(); anchor != anchors_.end(); ++anchor)
                {
                    if (isSecondAnchorFound
                        && std::abs(tgt::length(*anchor - grabbedAnchor_)) >= EPS
                        && std::abs(tgt::length(*anchor - grabbedSecondAnchor_)) >= EPS)
                    {
                        grabbedThirdAnchor_ = *anchor;
                        break;
                    }

                    if (!isSecondAnchorFound && std::abs(tgt::length(*anchor - grabbedAnchor_)) >= EPS) {
                        grabbedSecondAnchor_ = *anchor;
                        isSecondAnchorFound = true;
                    }
                }

                // Determine the cross product order of both direction vectors...
                tgt::vec3 planeNormal = tgt::normalize(planeNormal_.get());
                tgt::vec3 normal = tgt::normalize(tgt::cross(grabbedSecondAnchor_ - grabbedAnchor_, grabbedSecondAnchor_ - grabbedThirdAnchor_));
                if (std::abs(tgt::length(planeNormal - normal)) < EPS)
                //if (tgt::dot(planeNormal, normal) >= 0.f)
                    crossNewAnchorFirst_ = true;
                else
                    crossNewAnchorFirst_ = false;

                break;
            }

        // Check if a line was selected...
        for (std::vector<Line>::iterator it = lines_.begin(); it != lines_.end(); ++it)
            if (&(*it) == obj) {
                isLineGrabbed_ = true;
                grabbedLineBegin_ = *(it->first);
                grabbedLineEnd_ = *(it->second);

                tgt::vec3 begin = getWindowPos(grabbedLineBegin_);
                tgt::vec3 end = getWindowPos(grabbedLineEnd_);

                // t is almost correct for y-coordinate,
                // since the click event occured near the line center.
                float t = (e->coord().x - begin.x) / (end.x - begin.x);

                // Store modification line (with position vector given by click,
                // and direction vector equal to clipping plane normal vector)...
                grabbedLinePoint_ = grabbedLineBegin_ + t * (grabbedLineEnd_ - grabbedLineBegin_);
                grabbedLineDirection_ = tgt::normalize(planeNormal_.get());

                // Store clipping plane origin offset...
                grabbedPlaneOriginOffset_ = planePosition_.get();

                break;
            }

        // Is anchor or line grabbed?
        if (isAnchorGrabbed_ || isLineGrabbed_) {
            // Accept event, switch to interaction mode, and repaint...
            e->accept();
            QualityMode.requestQualityMode(VoreenQualityMode::RQ_INTERACTIVE, this);
            invalidate();
        }
    }

    // Mouse motion (applies modification in case an anchor or a line is grabbed)?
    if (e->action() & tgt::MouseEvent::MOTION) {
        // Is anchor or line grabbed?
        if (isAnchorGrabbed_ || isLineGrabbed_) {
            e->accept();

            // Is anchor grabbed?
            if (isAnchorGrabbed_) {
                tgt::vec3 mousePos = tgt::vec3(static_cast<float>(e->coord().x), e->viewport().y-static_cast<float>(e->coord().y), 0.0f);

                // Determine line vector equation for all 12 AABB edges...
                const int EDGES = 12;
                tgt::vec3 edgePoints[EDGES];
                tgt::vec3 edgeDirections[EDGES];

                tgt::Bounds bBox = inport_.getData()->getBoundingBox();
                tgt::vec3 llf = bBox.getLLF();
                tgt::vec3 urb = bBox.getURB();
                tgt::vec3 lrf = tgt::vec3(urb.x, llf.y, llf.z);
                tgt::vec3 lrb = tgt::vec3(urb.x, llf.y, urb.z);
                tgt::vec3 llb = tgt::vec3(llf.x, llf.y, urb.z);
                tgt::vec3 ulb = tgt::vec3(llf.x, urb.y, urb.z);
                tgt::vec3 ulf = tgt::vec3(llf.x, urb.y, llf.z);
                tgt::vec3 urf = tgt::vec3(urb.x, urb.y, llf.z);

                edgePoints[0] = edgePoints[1] = edgePoints[2] = llf;
                edgeDirections[0] = llb - llf;
                edgeDirections[1] = lrf - llf;
                edgeDirections[2] = ulf - llf;
                edgePoints[3] = edgePoints[4] = edgePoints[5] = urb;
                edgeDirections[3] = ulb - urb;
                edgeDirections[4] = urf - urb;
                edgeDirections[5] = lrb - urb;
                edgePoints[6] = edgePoints[7] = llb;
                edgeDirections[6] = ulb - llb;
                edgeDirections[7] = lrb - llb;
                edgePoints[8] = edgePoints[9] = ulf;
                edgeDirections[8] = ulb - ulf;
                edgeDirections[9] = urf - ulf;
                edgePoints[10] = edgePoints[11] = lrf;
                edgeDirections[10] = lrb - lrf;
                edgeDirections[11] = urf - lrf;

                // Determine edge on which the grabbed anchor lies...
                int edgeIndex = -1;
                float minDistance = 0.0f;
                for (int i = 0; i < EDGES; ++i) {
                    float distance = getPointLineDistance(grabbedAnchor_, edgePoints[i], edgeDirections[i]);
                    // No edge determined yet or distance lesser than minimum distance?
                    if (edgeIndex == -1 || distance < minDistance) {
                        edgeIndex = i;
                        minDistance = distance;
                    }
                }

                // Edge for grabbed anchor found?
                if (edgeIndex != -1) {
                    tgt::vec3 edgeStart = edgePoints[edgeIndex];
                    tgt::vec3 edgeEnd = edgePoints[edgeIndex] + edgeDirections[edgeIndex];
                    tgt::vec3 edgeDirection = edgeEnd - edgeStart;

                    tgt::vec3 edgeStartWP = getWindowPos(edgeStart);
                    tgt::vec3 edgeEndWP = getWindowPos(edgeEnd);
                    tgt::vec3 edgeDirectionWP = edgeEndWP - edgeStartWP;

                    // Determine position offset on the edge...
                    float t = getPointLineProjectionScalar(mousePos.xy(), edgeStartWP.xy(), edgeDirectionWP.xy());

                    // Is position offset near AABB corner according to JUMP_OVER_THRESHOLD?
                    const float JUMP_OVER_THRESHOLD = 0.01f;
                    if (t <= JUMP_OVER_THRESHOLD || t >= (1-JUMP_OVER_THRESHOLD)) {
                        // Determine the two additional edges which are connected with the nearest AABB corner..
                        const int CONNECTED_EDGES = 2;
                        tgt::vec3 connectedEdgesStart;
                        tgt::vec3 connectedEdgesEnd[CONNECTED_EDGES];

                        // Set appropriate edge start point of connected edges...
                        if (t >= (1-JUMP_OVER_THRESHOLD)) {
                            connectedEdgesStart = edgeEnd;
                        }
                        else if (t <= JUMP_OVER_THRESHOLD) {
                           connectedEdgesStart = edgeStart;
                        }

                        int connectedEdges = 0;
                        // Iterate over all AABB edges to find edge end points of connected edges...
                        for (int i = 0; i < EDGES; ++i) {
                            float edgeStartDistance = std::abs(tgt::length(edgePoints[i] - connectedEdgesStart));
                            float edgeEndDistance = std::abs(tgt::length((edgePoints[i] + edgeDirections[i]) - connectedEdgesStart));

                            // Is neither the AABB edge start point nor the AABB edge end point on the connected edge?
                            if (edgeStartDistance >= EPS && edgeEndDistance >= EPS)
                                continue;

                            // Is AABB edge start point equal to connected edge start point?
                            if (edgeStartDistance < EPS) {
                                connectedEdgesEnd[connectedEdges] = edgePoints[i] + edgeDirections[i];
                            }
                            // Is AABB edge end point equal to connected edge start point?
                            else if (edgeEndDistance < EPS)
                                connectedEdgesEnd[connectedEdges] = edgePoints[i];

                            // Are edges identical?
                            tgt::vec3 direction = connectedEdgesEnd[connectedEdges] - connectedEdgesStart;
                            if (std::abs(tgt::length(direction - edgeDirection)) < EPS
                                || std::abs(tgt::length(direction + edgeDirection)) < EPS)
                            {
                                continue;
                            }

                            // Were all connected edges found?
                            if (++connectedEdges == CONNECTED_EDGES)
                                    break;
                        }

                        // Jump over to another AABB edge, if mouse is nearer to this AABB edge than to the actual anchor AABB edge...
                        tgt::vec3 point = tgt::vec3(edgeStartWP.xy(), 0.0f);
                        tgt::vec3 direction = tgt::vec3(edgeDirectionWP.xy(), 0.0f);
                        float minDistance = getPointLineDistance(mousePos, point, direction);
                        bool tInRange = (0.0f <= t && t <= 1.0f ? true : false);
                        for (int i = 0; i < CONNECTED_EDGES; ++i) {
                            tgt::vec3 connectedEdgesDirection = connectedEdgesEnd[i] - connectedEdgesStart;
                            tgt::vec3 connectedEdgesStartWP = getWindowPos(connectedEdgesStart);
                            tgt::vec3 connectedEdgesEndWP = getWindowPos(connectedEdgesEnd[i]);
                            tgt::vec3 connectedEdgesDirectionWP = connectedEdgesEndWP - connectedEdgesStartWP;

                            point = tgt::vec3(connectedEdgesStartWP.xy(), 0.0f);
                            direction = tgt::vec3(connectedEdgesDirectionWP.xy(), 0.0f);
                            float distance = getPointLineDistance(mousePos, point, direction);
                            float edgeT = getPointLineProjectionScalar(mousePos.xy(), connectedEdgesStartWP.xy(), connectedEdgesDirectionWP.xy());

                            // Is scalar of the line vector equation of the actual anchor AABB edge out of range but
                            // the scalar of the line vector equation of the connected AABB edge in range
                            // or is the distance between mouse position and connected AABB edge lesser than
                            // actual minimum distance?
                            if ((!tInRange && 0.0f <= edgeT && edgeT <= 1.0f) || (distance < minDistance)) {
                                // Jump over to the connected AABB edge...
                                minDistance = distance;
                                t = edgeT;
                                tInRange = (0.0f <= t && t <= 1.0f ? true : false);

                                edgeStart = connectedEdgesStart;
                                edgeDirection = connectedEdgesDirection;
                            }
                        }
                    }

                    // Clamp t to [0, 1]...
                    t = std::min(std::max(t, 0.0f), 1.0f);

                    // Determine new position of grabbed anchor...
                    Anchor newAnchor = edgeStart + t * edgeDirection;

                    grabbedAnchor_ = newAnchor;

                    // Determine new clipping plane equation...
                    tgt::vec3 direction1 = grabbedSecondAnchor_ - newAnchor;
                    tgt::vec3 direction2 = grabbedSecondAnchor_ - grabbedThirdAnchor_;

                    tgt::vec3 planeNormal;
                    if (crossNewAnchorFirst_)
                        planeNormal = tgt::cross(direction1, direction2);
                    else
                        planeNormal = tgt::cross(direction2, direction1);

                    // Are plane direction vectors linearly independent?
                    if (std::abs(tgt::length(planeNormal)) >= EPS) {
                        // Manipulate clipping plane equation...
                        planeNormal = tgt::normalize(planeNormal);
                        float dist = tgt::dot(planeNormal, newAnchor - inport_.getData()->getBoundingBox().center());
                        planeNormal_.set(planeNormal);
                        planePosition_.set(dist);
                    }
                }
            }

            // Is a line grabbed?
            if (isLineGrabbed_) {
                // Determine mouse move offset...
                tgt::vec2 mousePos = tgt::vec2(static_cast<float>(e->coord().x), e->viewport().y-static_cast<float>(e->coord().y));
                tgt::vec3 begin = getWindowPos(grabbedLinePoint_);
                tgt::vec3 end = getWindowPos(grabbedLinePoint_ + grabbedLineDirection_);
                tgt::vec3 point = begin;
                tgt::vec3 direction = end - begin;
                float t = getPointLineProjectionScalar(mousePos, point.xy(), direction.xy());

                // Adapt clipping plane offset according to mouse move offset...
                tgt::vec4 plane = tgt::vec4(planeNormal_.get(), planePosition_.get());
                float length = tgt::length(plane.xyz());
                plane.w = grabbedPlaneOriginOffset_ + length * t;
                planeNormal_.set(plane.xyz());
                planePosition_.set(plane.w);
            }

            invalidate();
        }
    }

    // Mouse button released (indicates clipping plane manipulation end)?
    if (e->action() & tgt::MouseEvent::RELEASED) {
        if(isAnchorGrabbed_ || isLineGrabbed_) {
            // Accept mouse event, turn off interaction mode, ungrab anchor and lines, and repaint...
            e->accept();
            QualityMode.requestQualityMode(VoreenQualityMode::RQ_DEFAULT, this);
            isAnchorGrabbed_ = false;
            isLineGrabbed_ = false;
            invalidate();
        }
    }
}

void GeometryClippingWidget::renderPicking() {
    if (!idManager_)
        return;

    glPushAttrib(GL_ALL_ATTRIB_BITS);

    for (std::vector<Line>::iterator it = lines_.begin(); it != lines_.end(); ++it) {
        idManager_->setGLColor(&(*it));
        paintLine(*it);
    }

    for (std::vector<Anchor>::iterator it = anchors_.begin(); it != anchors_.end(); ++it) {
        idManager_->setGLColor(&(*it));
        paintAnchor(*it);
    }

    glPopAttrib();
}

void GeometryClippingWidget::render() {
    if (!renderGeometry_.get()) {
        return;
    }

    glPushAttrib(GL_ALL_ATTRIB_BITS);

    // Update anchors and lines, if plane changed...
    if (tgt::equal(actualPlane_, tgt::vec4(planeNormal_.get(), planePosition_.get())) != tgt::bvec4(true, true, true, true))
        updateAnchorsAndLines();

    // Render lines...
    glColor4fv(lineColor_.get().elem);
    for (std::vector<Line>::iterator it = lines_.begin(); it != lines_.end(); ++it)
        paintLine(*it);

    // Setup lighting...
    glEnable(GL_LIGHTING);
    glShadeModel(GL_SMOOTH);
    glLightfv(GL_LIGHT0, GL_POSITION, lightPosition_.get().elem);
    glEnable(GL_LIGHT0);

    tgt::vec3 spec(1.0f, 1.0f, 1.0f);
    glMaterialf( GL_FRONT_AND_BACK,    GL_SHININESS,    25.0f);
    glMaterialfv(GL_FRONT_AND_BACK,    GL_SPECULAR,    spec.elem);

    // Render plane on movement...
    glMaterialfv(GL_FRONT_AND_BACK,    GL_AMBIENT,        planeColor_.get().elem);
    glMaterialfv(GL_FRONT_AND_BACK,    GL_DIFFUSE,        planeColor_.get().elem);
    if (isAnchorGrabbed_ || isLineGrabbed_) {
        glBegin(GL_POLYGON);
        for (std::vector<Line>::iterator it = lines_.begin(); it != lines_.end(); ++it) {
            tgt::vec3 planeNormal = tgt::normalize(planeNormal_.get());
            glVertex3fv((*it->first + 0.01f * planeNormal).elem);
            glVertex3fv((*it->second + 0.01f * planeNormal).elem);
        }
        glEnd();
    }

    glMaterialfv(GL_FRONT_AND_BACK,    GL_AMBIENT,        tgt::vec4(1.0f, 0.0f, 0.0f, 1.0f).elem);
    glMaterialfv(GL_FRONT_AND_BACK,    GL_DIFFUSE,        tgt::vec4(1.0f, 0.0f, 0.0f, 1.0f).elem);

    // Render anchors...
    glMaterialfv(GL_FRONT_AND_BACK,    GL_AMBIENT,        anchorColor_.get().elem);
    glMaterialfv(GL_FRONT_AND_BACK,    GL_DIFFUSE,        anchorColor_.get().elem);
    for (std::vector<Anchor>::iterator it = anchors_.begin(); it != anchors_.end(); ++it) {
        // Is no anchor grabbed or the actual anchor not equal to the grabbed anchor?
        if (!isAnchorGrabbed_ || std::abs(tgt::length(*it - grabbedAnchor_)) >= EPS)
            paintAnchor(*it);
    }

    // Is anchor grabbed?
    if (isAnchorGrabbed_) {
        // Render grabbed anchor...
        glMaterialfv(GL_FRONT_AND_BACK,    GL_AMBIENT,        grabbedAnchorColor_.get().elem);
        glMaterialfv(GL_FRONT_AND_BACK,    GL_DIFFUSE,        grabbedAnchorColor_.get().elem);
        paintAnchor(grabbedAnchor_);
    }

    glPopAttrib();
}

void GeometryClippingWidget::paintAnchor(const Anchor& anchor) {
    GLUquadricObj* quadric = gluNewQuadric();
    MatStack.pushMatrix();
    tgt::vec3 planeNormal = tgt::normalize(planeNormal_.get());
    // Offset in clipping plane normal vector direction prevents z-fighting...
    tgt::vec3 position = anchor + 0.01f * planeNormal;
    MatStack.translate(position.x, position.y, position.z);
    gluSphere(quadric, 0.025f, 32, 32);
    MatStack.popMatrix();
    gluDeleteQuadric(quadric);
}

void GeometryClippingWidget::paintLine(const Line& line) {
    glLineWidth(width_.get());
    glBegin(GL_LINES);
    tgt::vec3 planeNormal = tgt::normalize(planeNormal_.get());
    // Offset in clipping plane normal vector direction prevents z-fighting...
    tgt::vec3 begin = *(line.first) + 0.01f * planeNormal;
    tgt::vec3 end = *(line.second) + 0.01f * planeNormal;
    glVertex3fv(begin.elem);
    glVertex3fv(end.elem);
    glEnd();
}

inline float GeometryClippingWidget::getPointLineDistance(const tgt::vec3& point, const tgt::vec3& linePositionVector, const tgt::vec3& lineDirectionVector) const {
    return tgt::length(tgt::cross(lineDirectionVector, point - linePositionVector)) / tgt::length(lineDirectionVector);
}

inline float GeometryClippingWidget::getPointLineProjectionScalar(const tgt::vec2& point, const tgt::vec2& linePositionVector, const tgt::vec2& lineDirectionVector) const {
    return tgt::dot(point - linePositionVector, lineDirectionVector) / tgt::lengthSq(lineDirectionVector);
}

bool GeometryClippingWidget::getIntersectionPointInAABB(tgt::vec3& intersectionPoint, const tgt::vec3& linePositionVector, const tgt::vec3& lineDirectionVector, const tgt::vec4& plane) {
    tgtAssert(inport_.getData(), "no input geometry");

    float denominator = tgt::dot(lineDirectionVector, plane.xyz());
    if (std::abs(denominator) < EPS)
        return false;

    float numerator = tgt::dot(linePositionVector, plane.xyz()) - plane.w;
    float t = -numerator / denominator;

    tgt::vec3 intersection = linePositionVector + t * lineDirectionVector;

    tgt::Bounds bBox = inport_.getData()->getBoundingBox();
    tgt::vec3 llf = bBox.getLLF();
    tgt::vec3 urb = bBox.getURB();

    // Is intersection point within AABB?
    if (llf.x <= intersection.x && intersection.x <= urb.x
        && llf.y <= intersection.y && intersection.y <= urb.y
        && llf.z <= intersection.z && intersection.z <= urb.z)
    {
        intersectionPoint = intersection;
        return true;
    }

    return false;
}

void GeometryClippingWidget::addAnchor(const tgt::vec3& linePositionVector, const tgt::vec3& lineDirectionVector, const tgt::vec4& plane) {
    tgt::vec3 intersectionPoint;
    // No intersection point withing the AABB?
    if (!getIntersectionPointInAABB(intersectionPoint, linePositionVector, lineDirectionVector, plane))
        return;

    // Omit anchor, if it already exist...
    for (std::vector<Anchor>::iterator it = anchors_.begin(); it != anchors_.end(); ++it)
        if (std::abs(tgt::length(*it - intersectionPoint)) < EPS)
            return;

    anchors_.push_back(intersectionPoint);
}

void GeometryClippingWidget::addLine(const tgt::vec4& planeAABB) {
    Anchor* begin = 0;
    Anchor* end = 0;

    for (std::vector<Anchor>::iterator it = anchors_.begin(); it != anchors_.end(); ++it) {
        // Is anchor on given plane?
        if (std::abs(tgt::dot(planeAABB.xyz(), *it) - planeAABB.w) < EPS) {
            if (!begin)
                begin = &(*it);
            else {
                end = &(*it);
                break;
            }
        }
    }

    // Was line start and end point found?
    if (begin && end)
        lines_.push_back(std::make_pair(begin, end));
}

void GeometryClippingWidget::updateAnchorsAndLines() {
    tgtAssert(inport_.getData(), "no input geometry");
    // Remove all old anchors and lines from ID manager and lists...
    if (idManager_) {
        for (std::vector<Anchor>::iterator it = anchors_.begin(); it != anchors_.end(); ++it)
            idManager_->deregisterObject(&(*it));
        for (std::vector<Line>::iterator it = lines_.begin(); it != lines_.end(); ++it)
            idManager_->deregisterObject(&(*it));
    }

    anchors_.clear();
    lines_.clear();

    tgt::Bounds bBox = inport_.getData()->getBoundingBox();
    tgt::vec3 llf = bBox.getLLF();
    tgt::vec3 urb = bBox.getURB();
    tgt::vec3 llb = tgt::vec3(llf.x, llf.y, urb.z);
    tgt::vec3 ulf = tgt::vec3(llf.x, urb.y, llf.z);
    tgt::vec3 lrf = tgt::vec3(urb.x, llf.y, llf.z);

    float length = tgt::length(planeNormal_.get());
    tgt::vec4 plane = tgt::vec4(planeNormal_.get() / length, planePosition_.get() + dot(planeNormal_.get(), inport_.getData()->getBoundingBox().center()));

    // Determine and add all anchors...
    addAnchor(llf, tgt::vec3( 1.0f,  0.0f,  0.0f), plane);
    addAnchor(llf, tgt::vec3( 0.0f,  1.0f,  0.0f), plane);
    addAnchor(llf, tgt::vec3( 0.0f,  0.0f,  1.0f), plane);
    addAnchor(urb, tgt::vec3(-1.0f,  0.0f,  0.0f), plane);
    addAnchor(urb, tgt::vec3( 0.0f, -1.0f,  0.0f), plane);
    addAnchor(urb, tgt::vec3( 0.0f,  0.0f, -1.0f), plane);
    addAnchor(llb, tgt::vec3( 1.0f,  0.0f,  0.0f), plane);
    addAnchor(llb, tgt::vec3( 0.0f,  1.0f,  0.0f), plane);
    addAnchor(ulf, tgt::vec3( 1.0f,  0.0f,  0.0f), plane);
    addAnchor(ulf, tgt::vec3( 0.0f,  0.0f,  1.0f), plane);
    addAnchor(lrf, tgt::vec3( 0.0f,  1.0f,  0.0f), plane);
    addAnchor(lrf, tgt::vec3( 0.0f,  0.0f,  1.0f), plane);

    // Less than 3 anchors (indicates that the clipping plane only touches the AABB)?
    if (anchors_.size() < 3)
        anchors_.clear();

    // Determine and add all lines according to connect the former determined anchors...
    addLine(tgt::vec4( 1.0f,  0.0f,  0.0f,  urb.x));
    addLine(tgt::vec4( 1.0f,  0.0f,  0.0f,  llf.x));
    addLine(tgt::vec4( 0.0f,  1.0f,  0.0f,  urb.y));
    addLine(tgt::vec4( 0.0f,  1.0f,  0.0f,  llf.y));
    addLine(tgt::vec4( 0.0f,  0.0f,  1.0f,  urb.z));
    addLine(tgt::vec4( 0.0f,  0.0f,  1.0f,  llf.z));

    // Reigster all new anchors and lines at the ID manager...
    if (idManager_) {
        for (std::vector<Anchor>::iterator it = anchors_.begin(); it != anchors_.end(); ++it)
            idManager_->registerObject(&(*it));
        for (std::vector<Line>::iterator it = lines_.begin(); it != lines_.end(); ++it)
            idManager_->registerObject(&(*it));
    }

    // Update actual plane equation to ensure that this update function is only called when necessary...
    actualPlane_ = tgt::vec4(planeNormal_.get(), planePosition_.get());

    // No lines found?
    if (lines_.size() == 0)
        return;

    //
    // Order lines and their anchors to CCW polygon vertex order which is expected by OpenGL...
    //

    // Sort lines to produce contiguous vertex order...
    bool reverseLastLine = false;
    for (std::vector<Line>::size_type i = 0; i < lines_.size() - 1; ++i) {
        for (std::vector<Line>::size_type j = i + 1; j < lines_.size(); ++j) {
            Anchor* connectionAnchor;
            if (reverseLastLine)
                connectionAnchor = lines_.at(i).first;
            else
                connectionAnchor = lines_.at(i).second;

            if (std::abs(tgt::length(*connectionAnchor - *(lines_.at(j).first))) < EPS) {
                std::swap(lines_.at(i + 1), lines_.at(j));
                reverseLastLine = false;
                break;
            }
            else if (std::abs(tgt::length(*connectionAnchor - *(lines_.at(j).second))) < EPS) {
                std::swap(lines_.at(i + 1), lines_.at(j));
                reverseLastLine = true;
                break;
            }
        }
    }

    std::vector<Anchor*> vertices;
    // Convert sorted line list to sorted vertex list...
    for (std::vector<Line>::size_type i = 0; i < lines_.size(); ++i) {
        bool reverseLine = i != 0 && std::abs(tgt::length(*vertices.at(vertices.size() - 1) - *(lines_.at(i).first))) >= EPS;

        Anchor* first = (reverseLine ? lines_.at(i).second : lines_.at(i).first);
        Anchor* second = (reverseLine ? lines_.at(i).first : lines_.at(i).second);

        if (i == 0)
            vertices.push_back(first);

        if (i < (vertices.size() - 1))
            vertices.push_back(second);
    }

    // Convert vertex order to counter clockwise if necessary...
    tgt::vec3 normal(0, 0, 0);
    for (std::vector<Anchor*>::size_type i = 0; i < vertices.size(); ++i)
        normal += tgt::cross(*vertices[i], *vertices[(i + 1) % vertices.size()]);
    normal = tgt::normalize(normal);

    if (tgt::dot(plane.xyz(), normal) < 0) {
        std::reverse(lines_.begin(), lines_.end());
        for (std::vector<Line>::iterator it = lines_.begin(); it != lines_.end(); ++it)
            std::swap(it->first, it->second);
    }
}

} //namespace voreen
