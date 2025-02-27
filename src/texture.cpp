#include "texture.h"
#include "color.h"

#include <assert.h>
#include <iostream>
#include <algorithm>

using namespace std;

namespace CS248 {

inline void uint8_to_float( float dst[4], unsigned char* src ) {
  uint8_t* src_uint8 = (uint8_t *)src;
  dst[0] = src_uint8[0] / 255.f;
  dst[1] = src_uint8[1] / 255.f;
  dst[2] = src_uint8[2] / 255.f;
  dst[3] = src_uint8[3] / 255.f;
}

inline void float_to_uint8( unsigned char* dst, float src[4] ) {
  uint8_t* dst_uint8 = (uint8_t *)dst;
  dst_uint8[0] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[0])));
  dst_uint8[1] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[1])));
  dst_uint8[2] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[2])));
  dst_uint8[3] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[3])));
}

void Sampler2DImp::generate_mips(Texture& tex, int startLevel) {

  // NOTE:
  // This starter code allocates the mip levels and generates a level
  // map by filling each level with placeholder data in the form of a
  // color that differs from its neighbours'. You should instead fill
  // with the correct data!

  // Advanced Task
  // Implement mipmap for trilinear filtering

  // check start level
  if ( startLevel >= tex.mipmap.size() ) {
    std::cerr << "Invalid start level";
  }

  // allocate sublevels
  int baseWidth  = tex.mipmap[startLevel].width;
  int baseHeight = tex.mipmap[startLevel].height;
  int numSubLevels = (int)(log2f( (float)max(baseWidth, baseHeight)));

  numSubLevels = min(numSubLevels, kMaxMipLevels - startLevel - 1);
  tex.mipmap.resize(startLevel + numSubLevels + 1);

  int width  = baseWidth;
  int height = baseHeight;
  for (int i = 1; i <= numSubLevels; i++) {

    MipLevel& level = tex.mipmap[startLevel + i];

    // handle odd size texture by rounding down
    width  = max( 1, width  / 2); assert(width  > 0);
    height = max( 1, height / 2); assert(height > 0);

    level.width = width;
    level.height = height;
    level.texels = vector<unsigned char>(4 * width * height);

  }

  // fill all 0 sub levels with interchanging colors
  Color colors[3] = { Color(1,0,0,1), Color(0,1,0,1), Color(0,0,1,1) };
  for(size_t i = 1; i < tex.mipmap.size(); ++i) {

    Color c = colors[i % 3];
    MipLevel& mip = tex.mipmap[i];

    for(size_t i = 0; i < 4 * mip.width * mip.height; i += 4) {
      float_to_uint8( &mip.texels[i], &c.r );
    }
  }

}

Color Sampler2DImp::sample_nearest(Texture& tex,
                                   float u, float v,
                                   int level) {

  // Task 4: Implement nearest neighbour interpolation

  // return magenta for invalid level
  if (level < 0 || level >= tex.mipmap.size()) {
    return Color(1,0,1,1);
  }

  MipLevel& mip = tex.mipmap[level];
  // Nearest neighbour: just round to the nearest pixel
  int x = static_cast<int>(u * mip.width);
  int y = static_cast<int>(v * mip.height);
  int idx = 4 * (x + y * mip.width);
  return Color(mip.texels[idx] / 255.0f,
               mip.texels[idx + 1] / 255.0f,
               mip.texels[idx + 2] / 255.0f,
               mip.texels[idx + 3] / 255.0f);
}

Color Sampler2DImp::sample_bilinear(Texture& tex,
                                    float u, float v,
                                    int level) {
    // return magenta for invalid level
    if (level < 0 || level >= tex.mipmap.size()) {
        return Color(1,0,1,1);
    }

    MipLevel& mip = tex.mipmap[level];

    // Convert to pixel coordinates and get fractional parts
    float texX = u * mip.width - 0.5f;
    float texY = v * mip.height - 0.5f;

    int x0 = static_cast<int>(floor(texX));
    int y0 = static_cast<int>(floor(texY));

    // Calculate fractional offsets
    float s = texX - x0;
    float t = texY - y0;

    // GL_CLAMP_TO_EDGE
    int width = static_cast<int>(mip.width);
    int height = static_cast<int>(mip.height);
    x0 = std::max(0, std::min(width - 1, x0));
    y0 = std::max(0, std::min(height - 1, y0));
    int x1 = std::max(0, std::min(width - 1, x0 + 1));
    int y1 = std::max(0, std::min(height - 1, y0 + 1));

    // Color of 4 nearest sample locations with texture values
    Color c00(&mip.texels[4 * (x0 + y0 * width)]);
    Color c10(&mip.texels[4 * (x1 + y0 * width)]);
    Color c01(&mip.texels[4 * (x0 + y1 * width)]);
    Color c11(&mip.texels[4 * (x1 + y1 * width)]);

    Color c0 = c00 * (1.0f - s) + c10 * s;
    Color c1 = c01 * (1.0f - s) + c11 * s;

    return c0 * (1.0f - t) + c1 * t;
}

Color Sampler2DImp::sample_trilinear(Texture& tex,
                                     float u, float v,
                                     float u_scale, float v_scale) {

  // Advanced Task
  // Implement trilinear filtering

  // return magenta for invalid level
  return Color(1,0,1,1);

}

} // namespace CS248
