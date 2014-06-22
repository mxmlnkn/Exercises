#pragma once

class Vec {
public:
    double x,y,z;
    Vec(void) : x(0),y(0),z(0) { }
    Vec(const Vec & v) : x(v.x),y(v.y),z(v.z) { }
    Vec(const double x, const double y, const double z) : x(x),y(y),z(z) { }
    double & operator[] (const uint16_t i);
    Vec& operator+= (const Vec & v);
    Vec operator+ (const Vec & v) const;
    Vec& operator-= (const Vec & v);
    Vec operator- (const Vec & v) const;
    double operator* (const Vec v) const;
    Vec& operator*= (const double a);
    Vec operator* (const double a) const;
    Vec& operator/= (const double a);
    Vec operator/ (const double a) const;
    bool operator== (const Vec v) const;
    bool operator!= (const Vec v) const;
    double norm() const;
    double norm2() const;
};
double & Vec::operator[] (const uint16_t i) {
    switch (i) {
        case 0: return x;
        case 1: return y;
        case 2: return z;
    }
    assert( i < 3 );
    return x;
}
double Vec::operator* (const Vec v) const {
    return this->x * v.x + this->y * v.y + this->z * v.z;
}

Vec& Vec::operator+= (const Vec & v) {
    this->x += v.x;
    this->y += v.y;
    this->z += v.z;
    return *this;
}
inline Vec Vec::operator+ (const Vec & v) const {
    Vec res = *this;
    res += v;
    return res;
}

inline Vec& Vec::operator-= (const Vec & v) {
    return (*this) += v * (-1);
}
inline Vec Vec::operator- (const Vec & v) const {
    return (*this)+( v*(-1) );
}

Vec& Vec::operator*= (const double a) {
    this->x *= a;
    this->y *= a;
    this->z *= a;
    return *this;
}
inline Vec Vec::operator* (const double a) const {
    Vec res = *this;
    res *= a;
    return res;
}

inline Vec& Vec::operator/= (const double a) {
    (*this) *= 1.0/a;
    return *this;
}
inline Vec Vec::operator/ (const double a) const {
    return (*this)*(1.0/a);
}

bool Vec::operator== (const Vec v) const {
    return ( (this->x == v.x) && (this->y == v.y) && (this->z == v.z) );
}
inline bool Vec::operator!= (const Vec v) const {
    return !((*this) == v);
}
double Vec::norm() const {
    double res = sqrt((*this).norm2());
    return res;
}
double Vec::norm2() const {
    assert( this->x == this->x );
    assert( this->y == this->y );
    assert( this->z == this->z );
    double res = this->x * this->x + this->y * this->y + this->z * this->z;
    assert( res == res );
    return res;
}
    

template<typename T>
Vec operator*(T const& scalar, Vec rhs)
{
    // scalar multiplication is commutative: s M = M s
    return rhs *= scalar; // calls rhs.operator*=(scalar);
}
