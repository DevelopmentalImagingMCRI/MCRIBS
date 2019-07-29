/*=========================================================================

  Library   : Image Registration Toolkit ()
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _COLOR_H
#define _COLOR_H

#define HAS_COLOR

#ifndef HAS_COLOR

typedef unsigned char Color;
typedef unsigned char ColorRGBA;

#else

class Color
{

public:

  unsigned char r;
  unsigned char g;
  unsigned char b;

  // Constructor (default)
  Color();

  // Constructor (copy)
  Color(const Color &);

  // Constructor (copy)
  Color(const ColorRGBA &);

  // Copy operators
  Color& operator= (const int &);

  // Copy operators
  Color& operator= (const Color &);

  // Copy operators
  Color& operator= (const ColorRGBA &);

  // Convert HSV color to RGB color
  void HSVtoRGB(double, double, double);

  // Set RGB color
  void RGBtoHSV(double, double, double);

};

inline Color::Color()
{
  r = 0;
  g = 0;
  b = 0;
}

inline Color::Color(const Color &c)
{
  r = c.r;
  g = c.g;
  b = c.b;
}

inline Color::Color(const ColorRGBA &c)
{
  r = c.r;
  g = c.g;
  b = c.b;
}

inline Color& Color::operator=(const Color &c)
{
  r = c.r;
  g = c.g;
  b = c.b;
  return *this;
}

inline Color& Color::operator=(const ColorRGBA &c)
{
  r = c.r;
  g = c.g;
  b = c.b;
  return *this;
}

inline Color& Color::operator=(const int &c)
{
  r = c;
  g = c;
  b = c;
  return *this;
}

#endif

#endif
