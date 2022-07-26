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

#include <mirtk/Image.h>

#include <mirtk/ColorRGBA.h>

void ColorRGBA::RGBtoHSV(double cr, double cg, double cb)
{
  r = static_cast<unsigned char>(round(255 * cr));
  g = static_cast<unsigned char>(round(255 * cg));
  b = static_cast<unsigned char>(round(255 * cb));
}

void ColorRGBA::HSVtoRGB(double h, double s, double v)
{
  int i;
  double f, p, q, t;

  if (s == 0) {
    if (h < 0) {
      r = static_cast<unsigned char>(round(255 * v));
      g = static_cast<unsigned char>(round(255 * v));
      b = static_cast<unsigned char>(round(255 * v));
    } else {
      std::cerr << "ColorRGBA::HSV: Undefined HSV color" << std::endl;
      exit(1);
    }
  } else {
    if (h == 1) h = 0;
    h = h * 6;
    i = (int)floor(h);
    f = h - i;
    p = v * (1 - s);
    q = v * (1 - (s * f));
    t = v * (1 - (s * (1 - f)));
    switch (i) {
    case 0:
      r = static_cast<unsigned char>(round(255 * v));
      g = static_cast<unsigned char>(round(255 * t));
      b = static_cast<unsigned char>(round(255 * p));
      break;
    case 1:
      r = static_cast<unsigned char>(round(255 * q));
      g = static_cast<unsigned char>(round(255 * v));
      b = static_cast<unsigned char>(round(255 * p));
      break;
    case 2:
      r = static_cast<unsigned char>(round(255 * p));
      g = static_cast<unsigned char>(round(255 * v));
      b = static_cast<unsigned char>(round(255 * t));
      break;
    case 3:
      r = static_cast<unsigned char>(round(255 * p));
      g = static_cast<unsigned char>(round(255 * q));
      b = static_cast<unsigned char>(round(255 * v));
      break;
    case 4:
      r = static_cast<unsigned char>(round(255 * t));
      g = static_cast<unsigned char>(round(255 * p));
      b = static_cast<unsigned char>(round(255 * v));
      break;
    case 5:
      r = static_cast<unsigned char>(round(255 * v));
      g = static_cast<unsigned char>(round(255 * p));
      b = static_cast<unsigned char>(round(255 * q));
      break;
    default:
      break;
    }
  }
}

