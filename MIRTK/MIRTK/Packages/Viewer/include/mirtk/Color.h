/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright (c) Imperial College London
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef _COLOR_H
#define _COLOR_H

#define HAS_COLOR
#ifndef HAS_COLOR

typedef unsigned char Color;
typedef unsigned char ColorRGBA;

#else

#include <mirtk/ViewerExport.h>


class MIRTK_Viewer_EXPORT Color
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
