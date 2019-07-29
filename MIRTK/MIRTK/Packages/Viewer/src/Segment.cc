/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <mirtk/Segment.h>


Segment::Segment()
{
  // Visibility
  _visible = false;

  // Label name
  _label = NULL;

  // Colour
  setColor(0, 0, 0);

  // Transperancy
  _trans = 0.0;
}

Segment::Segment(char* label, unsigned char cr, unsigned char cg, unsigned char cb, double t, int v)
{
  // Visibility
  _visible = v;

  // Label name
  _label = strdup(label);

  // Colour
  setColor(cr, cg, cb);

  // Transperancy
  _trans = t;
}

Segment::~Segment()
{
  if (_label != NULL) free(_label);
}

Segment& Segment::operator =(const Segment& s)
{
  strncpy(_hexColor, s._hexColor, HEX_LENGTH);
  _color = s._color;
  _label = strdup(s._label);
  _trans = s._trans;
  return *this;
}

void Segment::setLabel(char* label)
{
  if (_label != NULL) free(_label);
  if (label != NULL) {
    _label = strdup(label);
  } else {
    _label = NULL;
  }
}

void Segment::setColor(unsigned char r, unsigned char g, unsigned char b)
{
  _color.r = r;
  _color.g = g;
  _color.b = b;
  this->setHexColor();
}

void Segment::setHexColor(void)
{
  this->rgb2Hex(_color.r, _color.g, _color.b, _hexColor);
}

void Segment::setTrans(double t)
{
  _trans = t;
}

void Segment::setVisibility(int vis)
{
  _visible = vis;
}

void Segment::getColor(unsigned char *r, unsigned char *g, unsigned char *b) const
{
  *r = _color.r;
  *g = _color.g;
  *b = _color.b;
}

void Segment::getHex(char *hexColor) const
{
  strncpy(hexColor, _hexColor, HEX_LENGTH);
}

double Segment::getTrans() const
{
  return _trans;
}

char *Segment::getLabel() const
{
  return _label;
}

int Segment::getVisibility() const
{
  return _visible;
}


void Segment::rgb2Hex(int r, int g, int b, char* h)
{
  char* hr = int2Hex(r, 2);
  char* hg = int2Hex(g, 2);
  char* hb = int2Hex(b, 2);

  h[0] = hr[1];
  h[1] = hr[2];
  h[2] = hg[1];
  h[3] = hg[2];
  h[4] = hb[1];
  h[5] = hb[2];
  h[6] = '\0';

  delete [] hr;
  delete [] hg;
  delete [] hb;
}

char* Segment::int2Hex(int n, unsigned char round)
{
  unsigned char size = round;
  char i, *hex;
  const char *hexLookup = "0123456789ABCDEF";
  int temp = n;

  hex = new char[size+2];

  if (n<0) {
    hex[0]='-';
    n *= -1;
  } else {
    hex[0]=' ';
  }

  char mask = 0x000f;

  for (i=0; i<size; i++) {
    temp = n;
    temp >>=(4*i);
    temp &= mask;
    hex[size-i]=hexLookup[temp];
  }

  hex[size+1]= '\0';

  return hex;
}
