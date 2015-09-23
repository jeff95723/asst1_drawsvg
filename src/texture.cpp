#include "texture.h"

#include <assert.h>
#include <iostream>
#include <algorithm>

using namespace std;

namespace CMU462 {

    inline void uint8_to_float(float dst[4], unsigned char* src) {
        uint8_t* src_uint8 = (uint8_t *) src;
        dst[0] = src_uint8[0] / 255.f;
        dst[1] = src_uint8[1] / 255.f;
        dst[2] = src_uint8[2] / 255.f;
        dst[3] = src_uint8[3] / 255.f;
    }

    inline void float_to_uint8(unsigned char* dst, float src[4]) {
        uint8_t* dst_uint8 = (uint8_t *) dst;
        dst_uint8[0] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[0])));
        dst_uint8[1] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[1])));
        dst_uint8[2] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[2])));
        dst_uint8[3] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[3])));
    }

    void Sampler2DImp::generate_mips(Texture& tex, int startLevel) {

        // NOTE(sky):
        // The starter code allocates the mip levels and generates a level
        // map simply fills each level with a color that differs from its
        // neighbours'. The reference solution uses trilinear filtering
        // and it will only work when you have mipmaps.

        if (startLevel >= tex.mipmap.size()) {
            std::cerr << "Invalid start level";
        }

        // allocate sublevels
        int baseWidth = tex.mipmap[startLevel].width;
        int baseHeight = tex.mipmap[startLevel].height;
        int numSubLevels = (int) (log2f((float) max(baseWidth, baseHeight)));

        numSubLevels = min(numSubLevels, kMaxMipLevels - startLevel - 1);
        tex.mipmap.resize(startLevel + numSubLevels + 1);

        int width = baseWidth;
        int height = baseHeight;
        for (int i = 1; i <= numSubLevels; i++) {

            MipLevel& level = tex.mipmap[startLevel + i];

            // handle odd size texture by rounding down
            width = max(1, width / 2);
            assert(width > 0);
            height = max(1, height / 2);
            assert(height > 0);

            level.width = width;
            level.height = height;
            level.texels = vector<unsigned char>(4 * width * height);

        }

        // Use average colors from the previous level
        float buf[4];
        for (size_t i = 1; i < tex.mipmap.size(); ++i) {
            MipLevel& mip = tex.mipmap[i];
            MipLevel& prev_mip = tex.mipmap[i - 1];

            for (int x = 0; x < mip.width; x++) {
                for (int y = 0; y < mip.height; y++) {
                    Color tex_color = Color(0, 0, 0, 0);
                    for (int i = 0; i < 2; i++) {
                        for (int j = 0; j < 2; j++) {
                            int sx = 2 * x + i;
                            int sy = 2 * y + j;
                            int index = 4 * (sx + sy * prev_mip.width);
                            uint8_to_float(buf, &prev_mip.texels[index]);
                            tex_color += Color(buf[0] / 4, buf[1] / 4,
                                    buf[2] / 4, buf[3] / 4);
                        }
                    }

                    int index = 4 * (x + y * mip.width);
                    float_to_uint8(&mip.texels[index], &tex_color.r);
                }
            }
        }
    }

    Color Sampler2DImp::sample_nearest(Texture& tex, float u, float v,
            int level) {

        MipLevel &mip = tex.mipmap[level];
        int su = (int) round(u * mip.width);
        int sv = (int) round(v * mip.height);
        int index = 4 * (su + sv * mip.width);
        float color[4];
        uint8_to_float(color, &mip.texels[index]);
        return Color(color[0], color[1], color[2], color[3]);
    }

    Color Sampler2DImp::sample_bilinear(Texture& tex, float u, float v,
            int level) {

        MipLevel &mip = tex.mipmap[level];
        u = u * mip.width - 0.5;
        v = v * mip.height - 0.5;
        int index;
        float buf[4];
        Color color = Color(0, 0, 0, 0);

        int x = floor(u);
        int y = floor(v);
        float u_ratio = u - x;
        float v_ratio = v - y;
        float u_opposite = 1 - u_ratio;
        float v_opposite = 1 - v_ratio;

        /*  *|
         * -----
         *   |
         */
        index = 4 * (x + y * mip.width);
        uint8_to_float(buf, &mip.texels[index]);
        color += Color(buf[0], buf[1], buf[2], buf[3]) * u_opposite
                * v_opposite;

        /*   |*
         * -----
         *   |
         */
        index = 4 * (x + 1 + y * mip.width);
        uint8_to_float(buf, &mip.texels[index]);
        color += Color(buf[0], buf[1], buf[2], buf[3]) * u_ratio * v_opposite;

        /*   |
         * -----
         *  *|
         */
        index = 4 * (x + (y + 1) * mip.width);
        uint8_to_float(buf, &mip.texels[index]);
        color += Color(buf[0], buf[1], buf[2], buf[3]) * u_opposite * v_ratio;

        /*   |
         * -----
         *   |*
         */
        index = 4 * (x + 1 + (y + 1) * mip.width);
        uint8_to_float(buf, &mip.texels[index]);
        color += Color(buf[0], buf[1], buf[2], buf[3]) * u_ratio * v_ratio;

        return color;
    }

    Color Sampler2DImp::sample_trilinear(Texture& tex, float u, float v,
            float u_scale, float v_scale) {
        float level = log2f(sqrt(u_scale * u_scale + v_scale * v_scale));
        int lo_level = floor(level);
        int hi_level = ceil(level);
        float dlevel = level - lo_level;

        Color lo_color = sample_bilinear(tex, u, v);
        Color hi_color = sample_bilinear(tex, u, v);

        return lo_color * dlevel + hi_color * (1 - dlevel);
    }

} // namespace CMU462
