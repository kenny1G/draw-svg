#include "viewport.h"

#include "CS248.h"

namespace CS248 {

void ViewportImp::set_viewbox( float x, float y, float span ) {

  // Task 3 (part 2): 
  // Set svg to normalized device coordinate transformation. 
  // Your input arguments are defined as SVG canvas coordinates.
  this->x = x;
  this->y = y;
  this->span = span; 
    // Scale
  float scale_x = 1.0f / (2*span);
  float scale_y = scale_x;
    // Translate to align the top-left of the viewport to (0, 0)
  float trans_x = -(x - span);
  float trans_y = -(y - span);

  Matrix3x3 canvas_to_norm = Matrix3x3::identity();
  canvas_to_norm[0][0] = scale_x;                   
  canvas_to_norm[1][1] = scale_y;                  
  canvas_to_norm[2][0] = scale_x * trans_x;     
  canvas_to_norm[2][1] = scale_y * trans_y;     

  svg_2_norm = canvas_to_norm;
}

void ViewportImp::update_viewbox( float dx, float dy, float scale ) { 
  
  this->x -= dx;
  this->y -= dy;
  this->span *= scale;
  set_viewbox( x, y, span );
}

} // namespace CS248
