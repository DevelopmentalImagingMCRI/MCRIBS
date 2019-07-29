/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <mirtk/SegmentTable.h>


SegmentTable::SegmentTable()
{
}

SegmentTable::~SegmentTable()
{
}

void SegmentTable::Set(int id, char* label, unsigned char r, unsigned char g, unsigned char b, double trans, int vis)
{
  _entry[id].setLabel(label);
  _entry[id].setColor(r, g, b);
  _entry[id].setTrans(trans);
  _entry[id].setVisibility(vis);
}

void SegmentTable::SetLabel(int id, char* label)
{
  _entry[id].setLabel(label);
}

void SegmentTable::SetColor(int id, unsigned char red, unsigned char green, unsigned char blue)
{
  _entry[id].setColor(red, green, blue);
}

void SegmentTable::SetTrans(int id, double t)
{
  _entry[id].setTrans(t);
}

void SegmentTable::SetVisibility(int id, int vis)
{
  _entry[id].setVisibility(vis);
}

char *SegmentTable::Get(int id, unsigned char* r, unsigned char* g, unsigned char* b, double* trans, int* v) const
{
  // Get r,g,b
  _entry[id].getColor(r, g, b);

  // Get transparency
  *trans = _entry[id].getTrans();

  // Get visibility
  *v = _entry[id].getVisibility();

  return _entry[id].getLabel();
}

void SegmentTable::Clear()
{
  int i;

  for (i = 0; i < this->Size(); i++) {
    if (IsValid(i) == true) Clear(i);
  }
}

void SegmentTable::Clear(int id)
{
  _entry[id].setLabel(NULL);
  _entry[id].setColor(0, 0, 0);
  _entry[id].setTrans(0);
}

void SegmentTable::Read(char *name)
{
  double trans;
  char c, buffer[256];
  int i, n, r, g, b, id, vis;

  // Clear entries
  Clear();

  std::ifstream from(name);
  if (!from) {
    std::cerr << "SegmentTable::Read: Can't open file " << name << "\n";
    exit(1);
  }

  // Read keyword
  from >> buffer;
  if ((strcmp(buffer, "SegmentTable:") != 0) && (strcmp(buffer, "itkSegmentTable:") != 0)) {
    std::cerr << "SegmentTable::Read: Not a valid segment table" << std::endl;
    exit(1);
  }

  // Read number of valid entries
  from >> n;

  // Read entries
  for (i = 0; i < n; i++) {
    from >> id >> r >> g >> b >> trans >> vis;
    for (;;) {
      c = from.peek();
      if ((c == ' ') || (c == '\t')) {
        from.get(c);
      } else {
        break;
      }
    }
    from.getline(buffer, 255);
    _entry[id].setLabel(buffer);
    _entry[id].setColor(r, g, b);
    _entry[id].setTrans(trans);
    _entry[id].setVisibility(vis);
  }
}

void SegmentTable::Write(char *name)
{
  int id, n;
  unsigned char r, g, b;

  // Open file
  std::ofstream to(name);

  if (!to) {
    std::cerr << "SegmentTable::Write: Can't open file " << name << "\n";
    exit(1);
  }

  // Count number of valid entries
  n = 0;
  for (id = 0; id < this->Size(); id++) {
    if (IsValid(id) == true) n++;
  }

  // Write header
  to << "SegmentTable: " << n << std::endl;

  // Write entries
  for (id = 0; id < this->Size(); id++) {
    if (IsValid(id) == true) {
      _entry[id].getColor(&r, &g, &b);
      to << id << "\t" << int(r) << "\t" << int(g) << "\t" << int(b) << "\t" << _entry[id].getTrans()<< "\t" << _entry[id].getVisibility() << "\t" << _entry[id].getLabel() << std::endl;
    }
  }
}


