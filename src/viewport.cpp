#include "viewport.h"

#include "CMU462.h"

namespace CMU462 {

void ViewportImp::set_viewbox( float x, float y, float span ) {

  // Task 4 (part 2):
  // Set svg to normalized device coordinate transformation. Your input
  // arguments are defined as SVG canvans coordinates.
  this->x = x;
  this->y = y;
  this->span = span;
  float tlx = x - span,
        tly = y - span;
  Matrix3x3 transform_matrix = Matrix3x3::identity();
  transform_matrix(0,2) = - (tlx) / ( 2 * span);
  transform_matrix(1,2) = - (tly) / ( 2 * span);
  transform_matrix(0,0) = 1.0 / (2 * span);
  transform_matrix(1,1) = 1.0 / (2 * span);

  this->svg_2_norm = transform_matrix;
}

void ViewportImp::update_viewbox( float dx, float dy, float scale ) {

  this->x -= dx;
  this->y -= dy;
  this->span *= scale;
  set_viewbox( x, y, span );
}

} // namespace CMU462
