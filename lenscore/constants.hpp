#pragma once

namespace lens {

// M_PI is a POSIX extension, not standard C++; this project builds with
// CXX_EXTENSIONS OFF and targets MSVC too, so every header that needs pi
// pulls it from here instead of declaring its own copy. A single shared
// definition also avoids ODR collisions when two headers that each used to
// declare their own namespace-scope kPi are included together in one TU.
inline constexpr float kPi    = 3.14159265358979323846f;
inline constexpr float kTwoPi = 6.28318530717958647692f;

}  // namespace lens
