/*=========================================================================

  Library   : Image Registration Toolkit ()
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _LOOKUPTABLE_H
#define _LOOKUPTABLE_H

#include <mirtk/ColorRGBA.h>
#include <mirtk/Color.h>

typedef enum { ColorMode_Custom,
               ColorMode_HotMetal,
               ColorMode_Red,
               ColorMode_Green,
               ColorMode_Blue,
               ColorMode_Jacobian,
               ColorMode_JacobianExpansion,
               ColorMode_JacobianContraction,
               ColorMode_Luminance,
               ColorMode_InverseLuminance,
               ColorMode_Rainbow
             } ColorMode;

             
class LookupTable
{

	friend class RView;
	
  /// Min value of data
  int _minData;

  /// Max value of data
  int _maxData;

  /// Min value of display
  int _minDisplay;

  /// Max value of display
  int _maxDisplay;

  /// Color mode
  ColorMode _mode;

  /// Update lookup table
  void Update();

  /// Initialize lookup table with min and max values
  void Initialize(int, int);

  /// Set minimum and maximum display intensities
  void SetMinMaxDisplayIntensity(int,   int);

  /// Get minimum and maximum display intensities
  void GetMinMaxDisplayIntensity(int &, int &);

  /// Set minimum display intensity
  void SetMinDisplayIntensity(int);

  /// Set maximum display intensity
  void SetMaxDisplayIntensity(int);

  /// Get maximum display intensity
  int  GetMaxDisplayIntensity();

  /// Get minimum intensity
  int  GetMinIntensity();

  /// Get maximum intensity
  int  GetMaxIntensity();

 
public:

  /// Lookup table
  ColorRGBA *lookupTable;

  /// Constructor
  LookupTable(int = 0, int = 10000);

  /// Destructor
  ~LookupTable();

  /// Color scheme functions
  void SetColorModeToLuminance();
  void SetColorModeToInverseLuminance();
  void SetColorModeToRed();
  void SetColorModeToGreen();
  void SetColorModeToBlue();
  void SetColorModeToRainbow();
  void SetColorModeToHotMetal();
  void SetColorModeToJacobian();
  void SetColorModeToJacobianExpansion();
  void SetColorModeToJacobianContraction();

  /// Get color for given image value
  ColorRGBA &At(int);

  /// Get color for given image value
  ColorRGBA &operator ()(int);

  /// Get minimum display intensity
  int  GetMinDisplayIntensity();

  /// Return color scheme
  ColorMode GetColorMode();

  /// Read lookup table from file
  void Read(char *);

  /// Write lookup table to file
  void Write(char *);

};

inline ColorRGBA &LookupTable::At(int value)
{
  if      (value < _minData) return lookupTable[0];
  else if (value > _maxData) return lookupTable[_maxData];
  else                       return lookupTable[value];
}

inline ColorRGBA &LookupTable::operator ()(int value)
{
  return At(value);
}

inline int LookupTable::GetMinIntensity()
{
  return _minData;
}

inline int LookupTable::GetMaxIntensity()
{
  return _maxData;
}

inline int LookupTable::GetMinDisplayIntensity()
{
  return _minDisplay;
}

inline int LookupTable::GetMaxDisplayIntensity()
{
  return _maxDisplay;
}

inline void LookupTable::SetMinDisplayIntensity(int value)
{
  if (value > _maxData) {
    _minDisplay = _maxData;
  } else {
    _minDisplay = value;
  }
  this->Update();
}

inline void LookupTable::SetMaxDisplayIntensity(int value)
{
  if (value < _minData) {
    _maxDisplay = _minData;
  } else {
    _maxDisplay = value;
  }
  this->Update();
}

inline void LookupTable::SetMinMaxDisplayIntensity(int value1, int value2)
{
  if (value1 > _maxData) {
    _minDisplay = _maxData;
  } else {
    _minDisplay = value1;
  }
  this->Update();
  if (value2 < _minData) {
    _maxDisplay = _minData;
  } else {
    _maxDisplay = value2;
  }
  this->Update();
}

inline ColorMode LookupTable::GetColorMode()
{
  return _mode;
}


#endif
