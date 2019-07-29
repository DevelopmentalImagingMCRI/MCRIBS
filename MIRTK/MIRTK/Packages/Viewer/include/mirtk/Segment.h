/*=========================================================================

  Library   : Image Registration Toolkit ()
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _SEGMENT_H
#define _SEGMENT_H

#define HEX_LENGTH 7

#include <mirtk/Image.h>
#include <mirtk/ColorRGBA.h>
#include <mirtk/Color.h>

class Segment
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

  /// Transperancy
  double _trans;

  /// Visibility flag
  int _visible;

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
  void setVisibility(int);

  /// Return color
  void getColor(unsigned char*, unsigned char*, unsigned char*) const;

  /// Return color in hex
  void getHex(char *) const;

  /// Return transperancy
  double getTrans() const;

  /// Return label
  char *getLabel() const;

  /// Return if visible
  int   getVisibility() const;

  // General Methods
  void  rgb2Hex(int, int, int, char*);
  char* int2Hex(int, unsigned char);

};

#endif
