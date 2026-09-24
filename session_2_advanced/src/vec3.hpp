#pragma once

#include <cmath>

namespace axis {
inline constexpr int x = 0;
inline constexpr int y = 1;
inline constexpr int z = 2;
inline constexpr int count = 3;
}

struct Vec3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    double& operator[](int component) { return component == axis::x ? x : component == axis::y ? y : z; }
    double operator[](int component) const { return component == axis::x ? x : component == axis::y ? y : z; }
};

inline Vec3 operator+(const Vec3& a, const Vec3& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}

inline Vec3 operator-(const Vec3& a, const Vec3& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}

inline Vec3& operator+=(Vec3& a, const Vec3& b) {
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
    return a;
}

inline Vec3 operator*(double s, const Vec3& a) {
    return {s * a.x, s * a.y, s * a.z};
}

inline Vec3 operator*(const Vec3& a, double s) {
    return s * a;
}

inline Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

inline double dot(const Vec3& a, const Vec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline double norm(const Vec3& a) {
    return std::sqrt(dot(a, a));
}
