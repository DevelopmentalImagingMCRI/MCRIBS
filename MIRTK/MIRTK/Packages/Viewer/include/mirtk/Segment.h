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

#ifndef _SEGMENT_H
#define _SEGMENT_H

#define HEX_LENGTH 7

#include <mirtk/ViewerExport.h>
#include <mirtk/Image.h>
#include <mirtk/ColorRGBA.h>
#include <mirtk/Color.h>


class MIRTK_Viewer_EXPORT Segment
{

  friend class RView;
  friend class Viewer;

protected:

  /// Color
  Color _color;

  /// Color in hex format
  char _hexColor[HEX_LENGTH];

  /// Name of structure
  char *_label;

  /// Transparency
  double _trans;

  /// Visibility flag
  bool _visible;

  void setHexColor(void);

public:

  // Constructor (default)
  Segment();

  // Constructor (existing)
  Segment(char*, unsigned char, unsigned char, unsigned char, double, int = true);

  // Destructor
  virtual ~Segment(void);

  // Copy operator
  Segment& operator = (const Segment&);

  /// Set color
  void setColor(unsigned char, unsigned char, unsigned char);

  /// Set label
  void setLabel(char*);

  /// Set transperancy
  void setTrans(double);

  /// Set to visible
  void setVisibility(bool);

  /// Return color
  void getColor(unsigned char*, unsigned char*, unsigned char*) const;

  /// Return color in hex
  void getHex(char *) const;

  /// Return transperancy
  double getTrans() const;

  /// Return label
  char *getLabel() const;

  /// Return if visible
  bool getVisibility() const;

  // General Methods
  void  rgb2Hex(int, int, int, char*);
  char* int2Hex(int, unsigned char);

};

#endif
