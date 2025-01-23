#include "software_renderer.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "triangulation.h"

using namespace std;

namespace CS248 {

// Implements SoftwareRenderer //

// fill a sample location with color
void SoftwareRendererImp::fill_sample(int sx, int sy, const Color &color) {
  // Scale coordinates based on sample rate
  // For example, if sample_rate = 4 (2x2), each pixel has 4 samples
  int scaled_width = width * sample_rate;
  int scaled_height = height * sample_rate;

  // Check bounds against scaled dimensions
  if (sx < 0 || sx >= scaled_width)
    return;
  if (sy < 0 || sy >= scaled_height)
    return;

  Color sample_color;
  float inv255 = 1.0 / 255.0;
  sample_color.r = sexy_sample_buffer[4 * (sx + sy * scaled_width)] * inv255;
  sample_color.g =
      sexy_sample_buffer[4 * (sx + sy * scaled_width) + 1] * inv255;
  sample_color.b =
      sexy_sample_buffer[4 * (sx + sy * scaled_width) + 2] * inv255;
  sample_color.a =
      sexy_sample_buffer[4 * (sx + sy * scaled_width) + 3] * inv255;

  sample_color = alpha_blending(sample_color, color);

  sexy_sample_buffer[4 * (sx + sy * scaled_width)] =
      (uint8_t)(sample_color.r * 255);
  sexy_sample_buffer[4 * (sx + sy * scaled_width) + 1] =
      (uint8_t)(sample_color.g * 255);
  sexy_sample_buffer[4 * (sx + sy * scaled_width) + 2] =
      (uint8_t)(sample_color.b * 255);
  sexy_sample_buffer[4 * (sx + sy * scaled_width) + 3] =
      (uint8_t)(sample_color.a * 255);
}

// fill samples in the entire pixel specified by pixel coordinates
void SoftwareRendererImp::fill_pixel(int x, int y, const Color &color) {

  // Task 2: Re-implement this function

  // check bounds
  if (x < 0 || x >= width)
    return;
  if (y < 0 || y >= height)
    return;

  for (int sx = 0; sx < sample_rate; sx++) {
    for (int sy = 0; sy < sample_rate; sy++) {
      fill_sample(x * sample_rate + sx, y * sample_rate + sy, color);
    }
  }
}

void SoftwareRendererImp::draw_svg(SVG &svg) {

  // set top level transformation
  transformation = canvas_to_screen;

  // canvas outline
  Vector2D a = transform(Vector2D(0, 0));
  a.x--;
  a.y--;
  Vector2D b = transform(Vector2D(svg.width, 0));
  b.x++;
  b.y--;
  Vector2D c = transform(Vector2D(0, svg.height));
  c.x--;
  c.y++;
  Vector2D d = transform(Vector2D(svg.width, svg.height));
  d.x++;
  d.y++;

  svg_bbox_top_left = Vector2D(a.x + 1, a.y + 1);
  svg_bbox_bottom_right = Vector2D(d.x - 1, d.y - 1);

  // draw all elements
  for (size_t i = 0; i < svg.elements.size(); ++i) {
    draw_element(svg.elements[i]);
  }

  // draw canvas outline
  rasterize_line(a.x, a.y, b.x, b.y, Color::Black);
  rasterize_line(a.x, a.y, c.x, c.y, Color::Black);
  rasterize_line(d.x, d.y, b.x, b.y, Color::Black);
  rasterize_line(d.x, d.y, c.x, c.y, Color::Black);

  // resolve and send to pixel buffer
  resolve();
}

void SoftwareRendererImp::set_sample_rate(size_t sample_rate) {

  // Task 2:
  // You may want to modify this for supersampling support
  this->sample_rate = sample_rate;
  if (sexy_sample_buffer)
    delete[] sexy_sample_buffer;
  this->sexy_sample_buffer =
      new unsigned char[4 * width * height * sample_rate * sample_rate];
}

void SoftwareRendererImp::set_pixel_buffer(unsigned char *pixel_buffer,
                                           size_t width, size_t height) {
  this->pixel_buffer = pixel_buffer;
  this->width = width;
  this->height = height;
  if (sexy_sample_buffer)
    delete[] sexy_sample_buffer;
  this->sexy_sample_buffer =
      new unsigned char[4 * width * height * sample_rate * sample_rate];
}

void SoftwareRendererImp::draw_element(SVGElement *element) {

  // Task 3 (part 1):
  // Modify this to implement the transformation stack

  switch (element->type) {
  case POINT:
    draw_point(static_cast<Point &>(*element));
    break;
  case LINE:
    draw_line(static_cast<Line &>(*element));
    break;
  case POLYLINE:
    draw_polyline(static_cast<Polyline &>(*element));
    break;
  case RECT:
    draw_rect(static_cast<Rect &>(*element));
    break;
  case POLYGON:
    draw_polygon(static_cast<Polygon &>(*element));
    break;
  case ELLIPSE:
    draw_ellipse(static_cast<Ellipse &>(*element));
    break;
  case IMAGE:
    draw_image(static_cast<Image &>(*element));
    break;
  case GROUP:
    draw_group(static_cast<Group &>(*element));
    break;
  default:
    break;
  }
}

// Primitive Drawing //

void SoftwareRendererImp::draw_point(Point &point) {

  Vector2D p = transform(point.position);
  rasterize_point(p.x, p.y, point.style.fillColor);
}

void SoftwareRendererImp::draw_line(Line &line) {

  Vector2D p0 = transform(line.from);
  Vector2D p1 = transform(line.to);
  rasterize_line(p0.x, p0.y, p1.x, p1.y, line.style.strokeColor);
}

void SoftwareRendererImp::draw_polyline(Polyline &polyline) {

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

void SoftwareRendererImp::draw_rect(Rect &rect) {

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

void SoftwareRendererImp::draw_polygon(Polygon &polygon) {

  Color c;

  // draw fill
  c = polygon.style.fillColor;
  if (c.a != 0) {

    // triangulate
    vector<Vector2D> triangles;
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

void SoftwareRendererImp::draw_ellipse(Ellipse &ellipse) {

  // Advanced Task
  // Implement ellipse rasterization
}

void SoftwareRendererImp::draw_image(Image &image) {

  // Advanced Task
  // Render image element with rotation

  Vector2D p0 = transform(image.position);
  Vector2D p1 = transform(image.position + image.dimension);

  rasterize_image(p0.x, p0.y, p1.x, p1.y, image.tex);
}

void SoftwareRendererImp::draw_group(Group &group) {

  for (size_t i = 0; i < group.elements.size(); ++i) {
    draw_element(group.elements[i]);
  }
}

// Rasterization //

// The input arguments in the rasterization functions
// below are all defined in screen space coordinates

void SoftwareRendererImp::rasterize_point(float x, float y, Color color) {

  // fill in the nearest pixel
  int sx = (int)floor(x);
  int sy = (int)floor(y);

  // check bounds
  if (sx < 0 || sx >= width)
    return;
  if (sy < 0 || sy >= height)
    return;

  fill_pixel(sx, sy, color);
}

void SoftwareRendererImp::rasterize_line(float x0, float y0, float x1, float y1,
                                         Color color) {
  // High level: reduce 8 octant cases to only two (steep and shallow) by
  //    always drawing from right to left
  //
  // When calculating the error between the screen-coordinated line and the
  // conceptual 3d line, always assume the gradient is increasing so you only
  // need to use the general bresenham algorithm but keep track of the real
  // negative-ness of the line in y_step. Use this assumption to simplify and
  // generalize error calculation but use y_step when updating y when drawing
  // the line

  // Bresenham only liked integers
  int ix0 = (int)floor(x0);
  int iy0 = (int)floor(y0);
  int ix1 = (int)floor(x1);
  int iy1 = (int)floor(y1);

  // Handle straight lines
  if (ix0 == ix1) { // Vertical line
    int y_start = std::min(iy0, iy1);
    int y_end = std::max(iy0, iy1);
    for (int y = y_start; y <= y_end; y++) {
      rasterize_point(ix0, y, color);
    }
    return;
  }

  if (iy0 == iy1) { // Horizontal line
    int x_start = std::min(ix0, ix1);
    int x_end = std::max(ix0, ix1);
    for (int x = x_start; x <= x_end; x++) {
      rasterize_point(x, iy0, color);
    }
    return;
  }

  // Ensure we always draw left to right
  if (ix0 > ix1) {
    std::swap(ix0, ix1);
    std::swap(iy0, iy1);
  }

  int dx = ix1 - ix0;
  int dy = iy1 - iy0;
  int y_step = (dy > 0) ? 1 : -1;
  dy = abs(dy);

  int x = ix0;
  int y = iy0;

  // Handle different slope cases
  if (dx >= dy) {
    // Slope <= 1
    // More horizontal lines
    int eps = 2 * dy - dx;
    for (; x <= ix1; x++) {
      rasterize_point(x, y, color);
      // eps is tracking how far our screen space line thus far is from 3d space
      if (eps > 0) {
        y += y_step;
        eps += 2 * (dy - dx);
      } else {
        eps += 2 * dy;
      }
    }
  } else {
    // Slope > 1
    // More vertical lines
    int eps = 2 * dx - dy;
    int end_y = iy1;
    for (; y != end_y + y_step; y += y_step) {
      rasterize_point(x, y, color);
      if (eps > 0) {
        x++;
        eps += 2 * (dx - dy);
      } else {
        eps += 2 * dx;
      }
    }
  }

  // Advanced Task
  // Drawing Smooth Lines with Line Width
}

bool point_inside_triangle(float x0, float y0, float x1, float y1, float x2,
                           float y2, float sx, float sy) {
  // Calculate line equations L0, L1, L2 using V dot N
  // For each edge, calculate V (vector to test point) and N (2D normal to edge)

  // Edge 0: (x0,y0) to (x1,y1)
  float V0x = sx - x0, V0y = sy - y0;    // Vector to test point
  float N0x = y1 - y0, N0y = -(x1 - x0); // Normal to edge vector
  float L0 = V0x * N0x + V0y * N0y;      // V dot N

  // Edge 1: (x1,y1) to (x2,y2)
  float V1x = sx - x1, V1y = sy - y1;    // Vector to test point
  float N1x = y2 - y1, N1y = -(x2 - x1); // Normal to edge vector
  float L1 = V1x * N1x + V1y * N1y;      // V dot N

  // Edge 2: (x2,y2) to (x0,y0)
  float V2x = sx - x2, V2y = sy - y2;    // Vector to test point
  float N2x = y0 - y2, N2y = -(x0 - x2); // Normal to edge vector
  float L2 = V2x * N2x + V2y * N2y;      // V dot N

  // Point is inside if all line equations are negative
  return L0 < 0 && L1 < 0 && L2 < 0;
}

void SoftwareRendererImp::rasterize_triangle(float x0, float y0, float x1,
                                             float y1, float x2, float y2,
                                             Color color) {
  // Task 1:
  // Implement triangle rasterization
  // Draw triangle outline
  // bounding box coordinates
  float min_x = std::min({x0, x1, x2});
  float max_x = std::max({x0, x1, x2});
  float min_y = std::min({y0, y1, y2});
  float max_y = std::max({y0, y1, y2});

  for (int x = (int)min_x; x <= (int)max_x; x++) {
    for (int y = (int)min_y; y <= (int)max_y; y++) {
      // For each pixel, test all samples
      for (int sx = 0; sx < sample_rate; sx++) {
        for (int sy = 0; sy < sample_rate; sy++) {
          // Figure out where the sample is in pixel space
          float sample_x = x + (sx + 0.5f) / sample_rate;
          float sample_y = y + (sy + 0.5f) / sample_rate;

          if (point_inside_triangle(x0, y0, x1, y1, x2, y2, sample_x,
                                    sample_y)) {
            // Convert pixel coordinates to sample buffer coordinates
            int sample_buffer_x = x * sample_rate + sx;
            int sample_buffer_y = y * sample_rate + sy;
            fill_sample(sample_buffer_x, sample_buffer_y, color);
          }
        }
      }
    }
  }

  // Advanced Task
  // Implementing Triangle Edge Rules
}

void SoftwareRendererImp::rasterize_image(float x0, float y0, float x1,
                                          float y1, Texture &tex) {
  // Task 4:
  // Implement image rasterization
  // looping over x0 -> x1, y0 -> y1 cause of the specification:
  // "the image should cover all screen samples inside the specified rectangle"
  for (int x = static_cast<int>(std::floor(x0)); x <= static_cast<int>(std::floor(x1)); x++) {
    for (int y = static_cast<int>(std::floor(y0)); y <= static_cast<int>(std::floor(y1)); y++) {
      if (x < 0 || x >= width || y < 0 || y >= height) {
        continue;
      }

      float u = static_cast<float>(x) / static_cast<float>(width);
      float v = static_cast<float>(y) / static_cast<float>(height);
      Color color = sampler->sample_nearest(tex, u, v);
      fill_pixel(x, y, color);
    }
  }
}

// resolve samples to pixel buffer
void SoftwareRendererImp::resolve(void) {
  // For each pixel in the output buffer
  int scaled_width = width * sample_rate;
  for (int x = 0; x < width; x++) {
    for (int y = 0; y < height; y++) {
      Color avg_color(0,0,0,0);
      float inv255 = 1.0f / 255.0f;

      // Sum up all samples for this pixel
      for (int sx = 0; sx < sample_rate; sx++) {
        for (int sy = 0; sy < sample_rate; sy++) {
          int sample_x = x * sample_rate + sx;
          int sample_y = y * sample_rate + sy;
          int sample_idx = 4 * (sample_x + sample_y * scaled_width);

          Color sample_color(sexy_sample_buffer[sample_idx] * inv255,
                             sexy_sample_buffer[sample_idx + 1] * inv255,
                             sexy_sample_buffer[sample_idx + 2] * inv255,
                             sexy_sample_buffer[sample_idx + 3] * inv255);
          avg_color = avg_color + sample_color;
        }
      }

      float sample_count = sample_rate * sample_rate;
      avg_color = avg_color * (1.0f / sample_count);

      pixel_buffer[4 * (x + y * width)] = (uint8_t)(avg_color.r * 255);
      pixel_buffer[4 * (x + y * width) + 1] = (uint8_t)(avg_color.g * 255);
      pixel_buffer[4 * (x + y * width) + 2] = (uint8_t)(avg_color.b * 255);
      pixel_buffer[4 * (x + y * width) + 3] = (uint8_t)(avg_color.a * 255);
    }
  }
}

Color SoftwareRendererImp::alpha_blending(Color pixel_color, Color color) {
  float Er = color.r, Eg = color.g, Eb = color.b, Ea = color.a;
  float Cr = pixel_color.r, Cg = pixel_color.g, Cb = pixel_color.b, Ca = pixel_color.a;

  // src: https://www.w3.org/TR/SVGTiny12/painting.html#CompositingSimpleAlpha
  float Ca_prime = 1 - (1 - Ea) * (1 - Ca);
  float Cr_prime = (1 - Ea) * Cr + Er;
  float Cg_prime = (1 - Ea) * Cg + Eg;
  float Cb_prime = (1 - Ea) * Cb + Eb;

  return Color(Cr_prime, Cg_prime, Cb_prime, Ca_prime);
}

} // namespace CS248
