/*=========================================================================

  Library   : Image Registration Toolkit ()
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _COLORRGBA_H
#define _COLORRGBA_H

#define HAS_COLOR

#ifndef HAS_COLOR

typedef unsigned char ColorRGBA;

#else

class ColorRGBA
{

public:

  unsigned char r;
  unsigned char g;
  unsigned char b;
  float a;

  // Constructor (default)
  ColorRGBA();

  // Constructor (copy)
  ColorRGBA(const ColorRGBA &);

  // Copy operators
  ColorRGBA& operator= (const int &);

  // Copy operators
  ColorRGBA& operator= (const ColorRGBA &);

  // Convert HSV color to RGB color
  void HSVtoRGB(double, double, double);

  // Set RGB color
  void RGBtoHSV(double, double, double);

};

inline ColorRGBA::ColorRGBA()
{
  r = 0;
  g = 0;
  b = 0;
  a = 1;
}

inline ColorRGBA::ColorRGBA(const ColorRGBA &c)
{
  r = c.r;
  g = c.g;
  b = c.b;
  a = c.a;
}

inline ColorRGBA& ColorRGBA::operator=(const ColorRGBA &c)
{
  r = c.r;
  g = c.g;
  b = c.b;
  a = c.a;
  return *this;
}

inline ColorRGBA& ColorRGBA::operator=(const int &c)
{
  r = c;
  g = c;
  b = c;
  return *this;
}

#endif

#endif
