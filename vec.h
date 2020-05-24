#ifndef __vec__
#define __vec__

#include <cmath>
#include <iostream>
#include <cassert>

static const double pi = 4 * atan(1.0);

template<class T, int n> struct Vec;
template<class T, int n> T dot(const Vec<T,n>& u,const Vec<T,n>& v);

template<class T, int n>
struct Vec
{
    T x[n];

    Vec()
    {make_zero();}

    explicit Vec(const T& a)
    {assert(n == 1);x[0]=a;}

    Vec(const T& a, const T& b)
    {assert(n == 2);x[0]=a;x[1]=b;}

    Vec(const T& a, const T& b, const T& c)
    {assert(n == 3);x[0]=a;x[1]=b;x[2]=c;}

    template<class U>
    explicit Vec(const Vec<U,n>& v)
    {for(int i = 0; i < n; i++) x[i] = (T)v.x[i];}

    void make_zero()
    {fill(0);}

    void fill(T value)
    {for(int i = 0; i < n; i++) x[i] = value;}
    
    Vec& operator = (const Vec& v)
    {for(int i = 0; i < n; i++) x[i] = v.x[i]; return *this;}

    Vec& operator += (const Vec& v)
    {for(int i = 0; i < n; i++) x[i] += v.x[i]; return *this;}

    Vec& operator -= (const Vec& v)
    {for(int i = 0; i < n; i++) x[i] -= v.x[i]; return *this;}

    Vec& operator *= (const Vec& v)
    {for(int i = 0; i < n; i++) x[i] *= v.x[i]; return *this;}

    Vec& operator /= (const Vec& v)
    {for(int i = 0; i < n; i++) x[i] /= v.x[i]; return *this;}

    Vec& operator *= (const T& c)
    {for(int i = 0; i < n; i++) x[i] *= c; return *this;}

    Vec& operator /= (const T& c)
    {for(int i = 0; i < n; i++) x[i] /= c; return *this;}

    Vec operator + () const
    {return *this;}

    Vec operator - () const
    {Vec r; for(int i = 0; i < n; i++) r[i] = -x[i]; return r;}

    Vec operator + (const Vec& v) const
    {Vec r; for(int i = 0; i < n; i++) r[i] = x[i] + v.x[i]; return r;}

    Vec operator - (const Vec& v) const
    {Vec r; for(int i = 0; i < n; i++) r[i] = x[i] - v.x[i]; return r;}

    Vec operator * (const Vec& v) const
    {Vec r; for(int i = 0; i < n; i++) r[i] = x[i] * v.x[i]; return r;}

    Vec operator / (const Vec& v) const
    {Vec r; for(int i = 0; i < n; i++) r[i] = x[i] / v.x[i]; return r;}

    Vec operator * (const T& c) const
    {Vec r; for(int i = 0; i < n; i++) r[i] = x[i] * c; return r;}

    Vec operator / (const T& c) const
    {Vec r; for(int i = 0; i < n; i++) r[i] = x[i] / c; return r;}

    const T& operator[] (int i) const
    {return x[i];}

    T& operator[] (int i)
    {return x[i];}

    T magnitude_squared() const
    {return dot(*this, *this);}

    T magnitude() const
    {return sqrt(magnitude_squared());}

    // Be careful to handle the zero Vector gracefully
    Vec normalized() const
    {T mag = magnitude(); if(mag) return *this / mag; Vec r; r[0] = 1; return r;};
};

template <class T, int n>
inline Vec<T,n> operator * (const T& c, const Vec<T,n>& v)
{return v*c;}


template <class T, int n>
inline T dot(const Vec<T,n> & u, const Vec<T,n> & v)
{
    T r  =  0;
    for(int i = 0; i < n; i++) r += u.x[i] * v.x[i];
    return r;
}

template <class T >
inline Vec<T,3> cross(const Vec<T,3> & u, const Vec<T,3> & v)
{
    return Vec<T,3> (
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0]);
}

template<class T, int d>
inline Vec<T,d> componentwise_max(const Vec<T,d>& a, const Vec<T,d>& b)
{
    Vec<T,d> r;
    for(int i=0; i<d; i++) r[i] = std::max(a[i], b[i]);
    return r;
}

template<class T, int d>
inline Vec<T,d> componentwise_min(const Vec<T,d>& a, const Vec<T,d>& b)
{
    Vec<T,d> r;
    for(int i=0; i<d; i++) r[i] = std::min(a[i], b[i]);
    return r;
}

template <class T, int n>
std::ostream& operator << (std::ostream& out, const Vec<T,n> & u)
{
    for(int i = 0; i < n; i++)
    {
        if(i) out << ' ';
        out << u[i];
    }
    return out;
}

template <class T, int n>
std::istream& operator >> (std::istream& in, Vec<T,n> & u)
{
    for(int i = 0; i < n; i++)
    {
        in >> u[i];
    }
    return in;
}

typedef Vec<double,2> vec2;
typedef Vec<double,3> vec3;
typedef Vec<int,2> ivec2;
typedef Vec<int,3> ivec3;

#endif
