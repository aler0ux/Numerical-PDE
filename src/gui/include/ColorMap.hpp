#pragma once

#include <algorithm>

namespace pde {

struct RGB {
    unsigned char r, g, b;
};

// Viridis-like perceptual colormap for scalar fields. t is clamped to [0, 1].
inline RGB viridis(double t) {
    t = std::clamp(t, 0.0, 1.0);

    static const double cp[][3] = {
        {0.267, 0.005, 0.329}, {0.283, 0.141, 0.458}, {0.254, 0.265, 0.530},
        {0.207, 0.372, 0.553}, {0.164, 0.471, 0.558}, {0.128, 0.567, 0.551},
        {0.135, 0.659, 0.518}, {0.267, 0.749, 0.441}, {0.478, 0.821, 0.318},
        {0.741, 0.873, 0.150}, {0.993, 0.906, 0.144},
    };
    constexpr int n = 11;

    const double pos = t * (n - 1);
    int idx = static_cast<int>(pos);
    idx = std::clamp(idx, 0, n - 2);
    const double frac = pos - idx;

    auto lerp = [&](int c) { return cp[idx][c] + (cp[idx + 1][c] - cp[idx][c]) * frac; };

    return RGB{
        static_cast<unsigned char>(lerp(0) * 255.0),
        static_cast<unsigned char>(lerp(1) * 255.0),
        static_cast<unsigned char>(lerp(2) * 255.0),
    };
}

// Diverging blue-white-red colormap, useful for signed fields centered at 0.
inline RGB coolwarm(double t) {
    t = std::clamp(t, 0.0, 1.0);
    if (t < 0.5) {
        const double s = t / 0.5;
        return RGB{
            static_cast<unsigned char>((0.230 + s * (1.0 - 0.230)) * 255.0),
            static_cast<unsigned char>((0.299 + s * (1.0 - 0.299)) * 255.0),
            static_cast<unsigned char>((0.754 + s * (1.0 - 0.754)) * 255.0),
        };
    }
    const double s = (t - 0.5) / 0.5;
    return RGB{
        static_cast<unsigned char>((1.0 + s * (0.706 - 1.0)) * 255.0),
        static_cast<unsigned char>((1.0 + s * (0.016 - 1.0)) * 255.0),
        static_cast<unsigned char>((1.0 + s * (0.150 - 1.0)) * 255.0),
    };
}

} // namespace pde
