#include "software_renderer.h"

#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>

#include "triangulation.h"

using namespace std;

namespace CMU462 {

// Implements SoftwareRenderer //
    /* Takes two floats and swap the value of them */
    void swapVal(float &x, float &y);

    /* Checks if the given test point is in the triangle with
     * the endpoints of T0 T1 T2.
     */
    int inTriangle(Vector2D T0, Vector2D T1, Vector2D T2, Vector2D test);

    /* Checks if the test point is in the edge with the end points
     * T1 T2.
     */
    float ETest(Vector2D T1, Vector2D T2, Vector2D test);

    void SoftwareRendererImp::draw_svg(SVG& svg) {
        // Clear up target
        clear_target();

        // set top level transformation
        transformation = canvas_to_screen;

        // draw all elements
        for (size_t i = 0; i < svg.elements.size(); ++i) {
            draw_element(svg.elements[i]);
        }

        // draw canvas outline
        Vector2D a = transform(Vector2D(0, 0));
        a.x--;
        a.y++;
        Vector2D b = transform(Vector2D(svg.width, 0));
        b.x++;
        b.y++;
        Vector2D c = transform(Vector2D(0, svg.height));
        c.x--;
        c.y--;
        Vector2D d = transform(Vector2D(svg.width, svg.height));
        d.x++;
        d.y--;

        rasterize_line(a.x, a.y, b.x, b.y, Color::Black);
        rasterize_line(a.x, a.y, c.x, c.y, Color::Black);
        rasterize_line(d.x, d.y, b.x, b.y, Color::Black);
        rasterize_line(d.x, d.y, c.x, c.y, Color::Black);

        // resolve and send to render target
        resolve();

    }

    void SoftwareRendererImp::set_sample_rate(size_t sample_rate) {

        // Task 3:
        // You may want to modify this for supersampling support
        this->sample_rate = sample_rate;

        // Set up SSAA render target
        this->supersample_target_w = target_w * sample_rate;
        this->supersample_target_h = target_h * sample_rate;

        // Create super sample render target
        this->supersample_target = std::vector<unsigned char>(
                4 * supersample_target_w * supersample_target_h, 255);

    }

    void SoftwareRendererImp::set_render_target(unsigned char* render_target,
            size_t width, size_t height) {

        // Task 3:
        // You may want to modify this for supersampling support
        this->render_target = render_target;
        this->target_w = width;
        this->target_h = height;

        // Set up SSAA render target
        this->supersample_target_w = width * sample_rate;
        this->supersample_target_h = height * sample_rate;

        // Create super sample render target
        this->supersample_target = std::vector<unsigned char>(
                4 * supersample_target_w * supersample_target_h, 255);

    }

    void SoftwareRendererImp::draw_element(SVGElement* element) {

        // Save prev transformation
        Matrix3x3 prev_trans = transformation;
        // Update current transformation
        transformation = transformation * element->transform;

        switch (element->type) {
            case POINT:
                draw_point(static_cast<Point&>(*element));
                break;
            case LINE:
                draw_line(static_cast<Line&>(*element));
                break;
            case POLYLINE:
                draw_polyline(static_cast<Polyline&>(*element));
                break;
            case RECT:
                draw_rect(static_cast<Rect&>(*element));
                break;
            case POLYGON:
                draw_polygon(static_cast<Polygon&>(*element));
                break;
            case ELLIPSE:
                draw_ellipse(static_cast<Ellipse&>(*element));
                break;
            case IMAGE:
                draw_image(static_cast<Image&>(*element));
                break;
            case GROUP:
                draw_group(static_cast<Group&>(*element));
                break;
            default:
                break;
        }
        // restore previous transformation
        transformation = prev_trans;
    }

// Primitive Drawing //

    void SoftwareRendererImp::draw_point(Point& point) {

        Vector2D p = transform(point.position);
        rasterize_point(p.x, p.y, point.style.fillColor);

    }

    void SoftwareRendererImp::draw_line(Line& line) {

        Vector2D p0 = transform(line.from);
        Vector2D p1 = transform(line.to);
        rasterize_line(p0.x, p0.y, p1.x, p1.y, line.style.strokeColor);

    }

    void SoftwareRendererImp::draw_polyline(Polyline& polyline) {

        Color c = polyline.style.strokeColor;

        if (c.a != 0) {
            int nPoints = polyline.points.size();
            for (int i = 0; i < nPoints - 1; i++) {
                Vector2D p0 = transform(polyline.points[(i + 0) % nPoints]);
                Vector2D p1 = transform(polyline.points[(i + 1) % nPoints]);
                rasterize_line(p0.x, p0.y, p1.x, p1.y, c);
            }
        }
    }

    void SoftwareRendererImp::draw_rect(Rect& rect) {

        Color c;

        // draw as two triangles
        float x = rect.position.x;
        float y = rect.position.y;
        float w = rect.dimension.x;
        float h = rect.dimension.y;

        Vector2D p0 = transform(Vector2D(x, y));
        Vector2D p1 = transform(Vector2D(x + w, y));
        Vector2D p2 = transform(Vector2D(x, y + h));
        Vector2D p3 = transform(Vector2D(x + w, y + h));

        // draw fill
        c = rect.style.fillColor;
        if (c.a != 0) {
            rasterize_triangle(p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c);
            rasterize_triangle(p2.x, p2.y, p1.x, p1.y, p3.x, p3.y, c);
        }

        // draw outline
        c = rect.style.strokeColor;
        if (c.a != 0) {
            rasterize_line(p0.x, p0.y, p1.x, p1.y, c);
            rasterize_line(p1.x, p1.y, p3.x, p3.y, c);
            rasterize_line(p3.x, p3.y, p2.x, p2.y, c);
            rasterize_line(p2.x, p2.y, p0.x, p0.y, c);
        }

    }

    void SoftwareRendererImp::draw_polygon(Polygon& polygon) {

        Color c;

        // draw fill
        c = polygon.style.fillColor;
        if (c.a != 0) {

            // triangulate
            vector < Vector2D > triangles;
            triangulate(polygon, triangles);

            // draw as triangles
            for (size_t i = 0; i < triangles.size(); i += 3) {
                Vector2D p0 = transform(triangles[i + 0]);
                Vector2D p1 = transform(triangles[i + 1]);
                Vector2D p2 = transform(triangles[i + 2]);
                rasterize_triangle(p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c);
            }
        }

        // draw outline
        c = polygon.style.strokeColor;
        if (c.a != 0) {
            int nPoints = polygon.points.size();
            for (int i = 0; i < nPoints; i++) {
                Vector2D p0 = transform(polygon.points[(i + 0) % nPoints]);
                Vector2D p1 = transform(polygon.points[(i + 1) % nPoints]);
                rasterize_line(p0.x, p0.y, p1.x, p1.y, c);
            }
        }
    }

    void SoftwareRendererImp::draw_ellipse(Ellipse& ellipse) {

        // Extra credit

    }

    void SoftwareRendererImp::draw_image(Image& image) {

        Vector2D p0 = transform(image.position);
        Vector2D p1 = transform(image.position + image.dimension);

        rasterize_image(p0.x, p0.y, p1.x, p1.y, image.tex);
    }

    void SoftwareRendererImp::draw_group(Group& group) {

        for (size_t i = 0; i < group.elements.size(); ++i) {
            draw_element(group.elements[i]);
        }

    }

// Rasterization //

// The input arguments in the rasterization functions
// below are all defined in screen space coordinates

    void SoftwareRendererImp::rasterize_point(float x, float y, Color color) {

        // fill in the nearest pixel
        int sx = (int) floor(x);
        int sy = (int) floor(y);

        // check bounds
        if (sx < 0 || sx >= target_w)
            return;
        if (sy < 0 || sy >= target_h)
            return;

        int dx, dy, ssx, ssy;
        for (dx = 0; dx < sample_rate; dx++) {
            for (dy = 0; dy < sample_rate; dy++) {
                ssx = sample_rate * sx + dx;
                ssy = sample_rate * sy + dy;
                super_sample_rasterize_point(ssx, ssy, color);
            }
        }
    }

    /* Writes the point x, y in the super sample render target with
     * color color.
     */
    void SoftwareRendererImp::super_sample_rasterize_point(float x, float y,
            Color color) {
        float r = color.r, g = color.g, b = color.b, a = color.a;
        float canvas_r, canvas_g, canvas_b, canvas_a;
        float r_after, g_after, b_after, a_after;

        int sx = (int) floor(x), sy = (int) floor(y);

        int index = 4 * (sx + sy * supersample_target_w);

        // get canvas previous colors
        canvas_r = (supersample_target[index]) / 255.0;
        canvas_g = (supersample_target[index + 1]) / 255.0;
        canvas_b = (supersample_target[index + 2]) / 255.0;
        canvas_a = (supersample_target[index + 3]) / 255.0;

        // computing alpha composition
        a_after = a + canvas_a - a * canvas_a;
        if (a_after == 0) {
            r_after = 0;
            g_after = 0;
            b_after = 0;
        } else {
            r_after = (r * a + canvas_r * canvas_a * (1 - a)) / a_after;
            g_after = (g * a + canvas_g * canvas_a * (1 - a)) / a_after;
            b_after = (b * a + canvas_b * canvas_a * (1 - a)) / a_after;
        }

        supersample_target[index] = (uint8_t)(r_after * 255);
        supersample_target[index + 1] = (uint8_t)(g_after * 255);
        supersample_target[index + 2] = (uint8_t)(b_after * 255);
        supersample_target[index + 3] = (uint8_t)(a_after * 255);
    }

    void SoftwareRendererImp::rasterize_line(float x0, float y0, float x1,
            float y1, Color color) {
        float dx = x1 - x0, dy = y1 - y0;
        int bigSlope = abs(dy) > abs(dx);

        float vy;
        // Vertical Line
        if (x0 == x1) {
            for (vy = min(y0, y1); vy < max(y0, y1); vy++) {
                rasterize_point(x0, vy, color);
            }
            return;
        }

        if (bigSlope) {
            swapVal(x0, y0);
            swapVal(x1, y1);
        }

        if (x0 > x1) {
            swapVal(x0, x1);
            swapVal(y0, y1);
        }

        dx = x1 - x0, dy = y1 - y0;
        float x = x0, y = y0;
        float slope = dy / dx;
        float eps = 0;

        if (slope > 0) {
            for (x = x0; x <= x1; x++) {
                if (bigSlope) {
                    rasterize_point(y, x, color);
                } else {
                    rasterize_point(x, y, color);
                }
                if (eps + slope < 0.5) {
                    eps += slope;
                } else {
                    y++;
                    eps += (slope - 1);
                }
            }
        } else if (slope < 0) {
            for (x = x0; x <= x1; x++) {
                if (bigSlope) {
                    rasterize_point(y, x, color);
                } else {
                    rasterize_point(x, y, color);
                }
                if (eps + slope > -0.5) {
                    eps += slope;
                } else {
                    y--;
                    eps += (slope + 1);
                }
            }
        } else {
            // Horizontal line
            for (x = x0; x <= x1; x++) {
                rasterize_point(x, y0, color);
            }
        }
    }

    void SoftwareRendererImp::rasterize_triangle(float x0, float y0, float x1,
            float y1, float x2, float y2, Color color) {
        // Swap to counter clockwise
        Vector2D T0, T1, T2;
        if (x0 < x1) {
            T0 = Vector2D(x0, y0);
            T1 = Vector2D(x1, y1);
        } else {
            T0 = Vector2D(x1, y1);
            T1 = Vector2D(x0, y0);
        }

        T2 = Vector2D(x2, y2);
        if (ETest(T0, T1, T2) > 0) {
            T2 = T1;
            T1 = Vector2D(x2, y2);
        }

        float tlx, tly, brx, bry;
        float d = 1.0 / sample_rate;
        float sx, sy;
        int i, j;
        int x, y;
        tlx = min(min(x0, x1), x2);
        tly = min(min(y0, y1), y2);
        brx = max(max(x0, x1), x2);
        bry = max(max(y0, y1), y2);

        for (y = floor(tly); y < ceil(bry); y++) {
            for (x = floor(tlx); x < ceil(brx); x++) {
                // Go through all super sampling points
                for (i = 0; i < sample_rate; i++) {
                    for (j = 0; j < sample_rate; j++) {
                        sx = x + (i * d) + (d / 2);
                        sy = y + (j * d) + (d / 2);
                        if (inTriangle(T0, T1, T2, Vector2D(sx, sy))) {
                            super_sample_rasterize_point(sample_rate * sx,
                                    sample_rate * sy, color);
                        }
                    }
                }
            }
        }
    }

    void SoftwareRendererImp::rasterize_image(float x0, float y0, float x1,
            float y1, Texture& tex) {
        float x, y, u, v;
        float u_scale = 1.0 / transformation(0, 0);
        float v_scale = 1.0 / transformation(1, 1);
        Color color;
        for (x = floor(x0); x < ceil(x1); x++) {
            for (y = floor(y0); y < ceil(y1); y++) {
                // Get the corresponding point in texture coordinates
                u = (x + 0.5 - x0) / (x1 - x0);
                v = (y + 0.5 - y0) / (y1 - y0);
                color = this->sampler->sample_trilinear(tex, u, v, u_scale,
                        v_scale);

                rasterize_point(x + 0.5, y + 0.5, color);
            }
        }

    }

// resolve samples to render target
    void SoftwareRendererImp::resolve(void) {
        int x, y, sx, sy;
        int i, j, index;
        int divider = sample_rate * sample_rate;
        float r, g, b, a;

        for (x = 0; x < target_w; x++) {
            for (y = 0; y < target_h; y++) {
                // Re-init the color values
                r = 0;
                g = 0;
                b = 0;
                a = 0;
                for (i = 0; i < sample_rate; i++) {
                    for (j = 0; j < sample_rate; j++) {
                        sx = sample_rate * x + i;
                        sy = sample_rate * y + j;
                        index = 4 * (sx + sy * supersample_target_w);
                        r += supersample_target[index] / 255.0;
                        g += supersample_target[index + 1] / 255.0;
                        b += supersample_target[index + 2] / 255.0;
                        a += supersample_target[index + 3] / 255.0;
                    }
                }

                // Take the average and write to the render target
                index = 4 * (x + y * target_w);
                render_target[index] = (uint8_t)((r / divider) * 255);
                render_target[index + 1] = (uint8_t)((g / divider) * 255);
                render_target[index + 2] = (uint8_t)((b / divider) * 255);
                render_target[index + 3] = (uint8_t)((a / divider) * 255);
            }
        }

        return;

    }

    void swapVal(float &x, float &y) {
        float tmp = y;
        y = x;
        x = tmp;
    }

    float ETest(Vector2D T1, Vector2D T2, Vector2D test) {
        Vector2D dP = T2 - T1;

        Vector2D dTest = T1 - test;   // Offset vector from test pointing to T1.

        // Computes the 2D perpendicular direction.
        Vector2D dP_perp = Vector2D(dP.y, -dP.x);

        return dot(dTest, dP_perp);
    }

    int inTriangle(Vector2D T0, Vector2D T1, Vector2D T2, Vector2D test) {
        return ((ETest(T0, T1, test) <= 0) && (ETest(T1, T2, test) <= 0)
                && (ETest(T2, T0, test) <= 0));
    }

} // namespace CMU462
