/*
 * spectrum.h
 * dnelson
 */
 
#ifndef AREPO_RT_SPECTRUM_H
#define AREPO_RT_SPECTRUM_H

#include "ArepoRT.h"

inline void XYZToRGB(const float xyz[3], float rgb[3]) {
  rgb[0] =  3.240479f*xyz[0] - 1.537150f*xyz[1] - 0.498535f*xyz[2];
  rgb[1] = -0.969256f*xyz[0] + 1.875991f*xyz[1] + 0.041556f*xyz[2];
  rgb[2] =  0.055648f*xyz[0] - 0.204043f*xyz[1] + 1.057311f*xyz[2];
}
inline void RGBToXYZ(const float rgb[3], float xyz[3]) {
  xyz[0] = 0.412453f*rgb[0] + 0.357580f*rgb[1] + 0.180423f*rgb[2];
  xyz[1] = 0.212671f*rgb[0] + 0.715160f*rgb[1] + 0.072169f*rgb[2];
  xyz[2] = 0.019334f*rgb[0] + 0.119193f*rgb[1] + 0.950227f*rgb[2];
}

template <int nSamples> class CoefficientSpectrum {
public:
  // construction
  CoefficientSpectrum(float v = 0.0f) {
    for (int i = 0; i < nSamples; ++i)
      c[i] = v;
  }

  // operators
  CoefficientSpectrum &operator+=(const CoefficientSpectrum &s2) {
    for (int i = 0; i < nSamples; ++i)
      c[i] += s2.c[i];
    return *this;
  }
  CoefficientSpectrum operator+(const CoefficientSpectrum &s2) const {
    CoefficientSpectrum ret = *this;
    for (int i = 0; i < nSamples; ++i)
      ret.c[i] += s2.c[i];
    return ret;
  }
  CoefficientSpectrum operator-(const CoefficientSpectrum &s2) const {
    CoefficientSpectrum ret = *this;
    for (int i = 0; i < nSamples; ++i)
      ret.c[i] -= s2.c[i];
    return ret;
  }
  CoefficientSpectrum operator/(const CoefficientSpectrum &s2) const {
    CoefficientSpectrum ret = *this;
    for (int i = 0; i < nSamples; ++i)
      ret.c[i] /= s2.c[i];
    return ret;
  }
  CoefficientSpectrum operator*(const CoefficientSpectrum &sp) const {
    CoefficientSpectrum ret = *this;
    for (int i = 0; i < nSamples; ++i)
      ret.c[i] *= sp.c[i];
    return ret;
  }
  CoefficientSpectrum &operator*=(const CoefficientSpectrum &sp) {
    for (int i = 0; i < nSamples; ++i)
      c[i] *= sp.c[i];
    return *this;
  }
  CoefficientSpectrum operator*(float a) const {
    CoefficientSpectrum ret = *this;
    for (int i = 0; i < nSamples; ++i)
      ret.c[i] *= a;
    return ret;
  }
  CoefficientSpectrum &operator*=(float a) {
    for (int i = 0; i < nSamples; ++i)
      c[i] *= a;
    return *this;
  }
  friend inline
  CoefficientSpectrum operator*(float a, const CoefficientSpectrum &s) {
    return s * a;
  }
  CoefficientSpectrum operator/(float a) const {
    CoefficientSpectrum ret = *this;
    for (int i = 0; i < nSamples; ++i)
      ret.c[i] /= a;
    return ret;
  }
  CoefficientSpectrum &operator/=(float a) {
    for (int i = 0; i < nSamples; ++i)
      c[i] /= a;
    return *this;
  }
  bool operator==(float a) const {
    for (int i = 0; i < nSamples; ++i)
      if (c[i] != a) return false;
    return true;
  }
  bool operator==(const CoefficientSpectrum &sp) const {
    for (int i = 0; i < nSamples; ++i)
      if (c[i] != sp.c[i]) return false;
    return true;
  }
  bool operator!=(const CoefficientSpectrum &sp) const {
    return !(*this == sp);
  }
  
  CoefficientSpectrum operator-() const {
    CoefficientSpectrum ret;
    for (int i = 0; i < nSamples; ++i)
      ret.c[i] = -c[i];
    return ret;
  }
  friend CoefficientSpectrum Exp(const CoefficientSpectrum &s) {
    CoefficientSpectrum ret;
    for (int i = 0; i < nSamples; ++i)
      ret.c[i] = expf(s.c[i]);
    return ret;
  }
  
  // data
protected:
  float c[nSamples];
};

class RGBSpectrum : public CoefficientSpectrum<3> {
  using CoefficientSpectrum<3>::c;
public:
  // construction
  RGBSpectrum(float v = 0.0f) : CoefficientSpectrum<3>(v) {
  }
  RGBSpectrum(const CoefficientSpectrum<3> &v)
      : CoefficientSpectrum<3>(v) {
  }
  RGBSpectrum(const RGBSpectrum &s) {
    *this = s;
  }
  
  // conversions
  static RGBSpectrum FromRGB(const float rgb[3]) {
    RGBSpectrum s;
    s.c[0] = rgb[0];
    s.c[1] = rgb[1];
    s.c[2] = rgb[2];
    return s;
  }
  void ToRGB(float *rgb) const {
    rgb[0] = c[0];
    rgb[1] = c[1];
    rgb[2] = c[2];
  }
  void ToXYZ(float xyz[3]) const {
    RGBToXYZ(c, xyz);
  }
  const RGBSpectrum &ToRGBSpectrum() const {
    return *this;
  }
  
  // TODO: lookup from table
  static RGBSpectrum FromNamed(const string &name) {
    RGBSpectrum s;
    if (name == "red") {
      s.c[0] = 0.1f;
      s.c[1] = 0.0f;
      s.c[2] = 0.0f;
    } else if (name == "green") {
      s.c[0] = 0.0f;
      s.c[1] = 0.1f;
      s.c[2] = 0.0f;
    } else if (name == "blue") {
      s.c[0] = 0.0f;
      s.c[1] = 0.0f;
      s.c[2] = 0.1f;
    } else if (name == "black") {
      s.c[0] = 0.0f;
      s.c[1] = 0.0f;
      s.c[2] = 0.0f;
    } else if (name == "white") {
      s.c[0] = 0.1f;
      s.c[1] = 0.1f;
      s.c[2] = 0.1f;
    }
    return s;
  }
  
  // components
  float r() { return c[0]; }
  float g() { return c[1]; }
  float b() { return c[2]; }
  
  float y() const {
    const float YWeight[3] = { 0.212671f, 0.715160f, 0.072169f };
    return YWeight[0] * c[0] + YWeight[1] * c[1] + YWeight[2] * c[2];
  }

};
  
#endif
