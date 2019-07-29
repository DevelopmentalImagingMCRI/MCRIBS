/*=========================================================================

  Library   : Image Registration Toolkit ()
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _RVIEW_H
#define _RVIEW_H

#include <iostream>

#include <mirtk/IOConfig.h>

typedef enum { View_A,
               View_B,
               View_Checkerboard,
               View_Subtraction,
               View_HShutter,
               View_VShutter,
               View_AoverB,
               View_BoverA
             } RViewMode;

typedef enum { NoneDef,
               Displacement,
               Jacobian,
               Jacobian_Expansion,
               Jacobian_Contraction
             } DeformationProperty;

typedef enum { Viewer_XY, Viewer_XZ, Viewer_YZ, Viewer_None } ViewerMode;

typedef enum { CrossHair, CursorX, CursorV, CursorBar } CursorMode;

typedef enum { Neurological, Radiological, Native } DisplayMode;

typedef enum { RegionGrowing2D, RegionGrowing3D } RegionGrowingMode;

#define DEFORMATION_DISPLACEMENT_MIN 0
#define DEFORMATION_DISPLACEMENT_MAX 100
#define DEFORMATION_JACOBIAN_MIN 10
#define DEFORMATION_JACOBIAN_MAX 400

// Function keys
#define KEY_F1                     1
#define KEY_F2                     2
#define KEY_F3                     3
#define KEY_F4                     4
#define KEY_F5                     5
#define KEY_F6                     6
#define KEY_F7                     7
#define KEY_F8                     8
#define KEY_F9                     9
#define KEY_F10                    10
#define KEY_F11                    11
#define KEY_F12                    12

// Directional keys
#define KEY_LEFT                   100
#define KEY_UP                     101
#define KEY_RIGHT                  102
#define KEY_DOWN                   103
#define KEY_PAGE_UP                104

// Mouse buttons.
#define MOUSE_LEFT_BUTTON                0
#define MOUSE_MIDDLE_BUTTON              1
#define MOUSE_RIGHT_BUTTON               2

// Mouse button  state.
#define MOUSE_DOWN                       0
#define MOUSE_UP                         1

#ifndef IMPERIAL

//#include <mirtk/TransformationCollection.h>

#define AFFDTransformation "AdaptiveFreeFormTransformation"
#define MFFDTransformation "mirtk::TransformationCollection"

#define BSplineFreeFormTransformation AdaptiveFreeFormTransformation
#define MFreeFormTransformation mirtk::TransformationCollection
#define MFreeFormRegistration   AdaptiveFreeFormRegistration
#define MFreeFormRegistration2D AdaptiveFreeFormRegistration2D

#endif

#include <list>

#if MIRTK_IO_WITH_VTK && defined(HAVE_VTK)
#include <vtkPointSet.h>
#endif

#include <mirtk/ImageTransformation.h>
#include <mirtk/MultiLevelFreeFormTransformation.h>

#define MAX_SEGMENTS 256

#define MAX_NUMBER_OF_OBJECTS 40

#include <mirtk/SegmentTable.h>

#include <mirtk/LookupTable.h>
#include <mirtk/Viewer.h>
#include <mirtk/RViewConfig.h>
#include <mirtk/HistogramWindow.h>
#include <mirtk/VoxelContour.h>

class VoxelContour;

class RView
{

protected:

  /// Friends
  friend class Viewer;
  friend class LookupTable;
  friend class VoxelContour;
  friend class SegmentationEditor;

  /// Number of image viewers
  int _NoOfViewers;

  /// Image viewer for target image
  Viewer **_viewer;

  /// Whether the given viewer displays the source image instead
  bool *_isSourceViewer;

  /// Target image
  mirtk::Image *_targetImage;

  /// Source image
  mirtk::Image *_sourceImage;

  /// Segmentation image
  mirtk::GreyImage *_segmentationImage;

  /// Segment Table
  SegmentTable *_segmentTable;

  /// Color Table
  Color *_segmentColorTable;

  /// Transformation for reslicing of target image
  mirtk::Transformation *_targetTransform;

  /// Transformation for reslicing of source image
  mirtk::Transformation *_sourceTransform;

  /// Transformation for reslicing of segmentation image
  mirtk::Transformation *_segmentationTransform;

  /// Transformation for reslicing of selection image
  mirtk::Transformation *_selectionTransform;

  /// Transformation filter for reslicing of target image
  mirtk::ImageTransformation **_targetTransformFilter;

  /// Transformation filter for reslicing of source image
  mirtk::ImageTransformation **_sourceTransformFilter;

  /// Displacement field used to cache source transformation
  mirtk::ImageTransformationCache _sourceTransformCache;

  /// Whether to cache displacements or not
  int _CacheDisplacements;

  /// Transformation filter for reslicing of segmentation image
  mirtk::ImageTransformation **_segmentationTransformFilter;

  /// Transformation filter for reslicing of selection image
  mirtk::ImageTransformation **_selectionTransformFilter;

  /// Target image
  mirtk::GreyImage **_targetImageOutput;

  /// Source image
  mirtk::GreyImage **_sourceImageOutput;

  /// Segmentation image
  mirtk::GreyImage **_segmentationImageOutput;

  /// Selection image
  mirtk::GreyImage **_selectionImageOutput;

  /// Target landmarks (pointset)
  mirtk::PointSet _targetLandmarks;

  /// Source landmarks (pointset)
  mirtk::PointSet _sourceLandmarks;

  /// Contour
  VoxelContour _voxelContour;

  /// Contour viewer
  int _contourViewer;

  /// Contour viewer mode
  ViewerMode _contourViewerMode;

#if MIRTK_IO_WITH_VTK && defined(HAVE_VTK)
  /// Number of vtk objects
  int _NoOfObjects;

  /// Object (surface or mesh)
  vtkPointSet *_Object[MAX_NUMBER_OF_OBJECTS];

  /// Flag for display object as a movie
  int _ObjectMovie;

#endif

  /// Combined source and target images in OpenGL format
  Color **_drawable;

  /// Color lookup table for target image
  LookupTable *_targetLookupTable;

  /// Color lookup table for source image
  LookupTable *_sourceLookupTable;

  /// Color lookup table for subtraction of target and source image
  LookupTable *_subtractionLookupTable;

  /// Target value range
  double _targetMin, _targetMax;

  /// Source value range
  double _sourceMin, _sourceMax;

  /// Subtraction value range
  double _subtractionMin, _subtractionMax;

  /// Target display value range
  double _targetDisplayMin, _targetDisplayMax;
  
  /// Source display value range
  double _sourceDisplayMin, _sourceDisplayMax;

  /// Subtraction display value range
  double _subtractionDisplayMin, _subtractionDisplayMax;

  /// Target frame
  int _targetFrame;

  /// Source frame
  int _sourceFrame;

  /// Flag to indicate whether target image must be updated
  int _targetUpdate;

  /// Flag to indicate whether source image must be updated
  int _sourceUpdate;

  /// Flag to indicate whether source image must be updated
  int _segmentationUpdate;

  /// Flag to indicate whether selection image must be updated
  int _selectionUpdate;

  /// Width of viewer  (in pixels)
  int _screenX;

  /// Height of viewer (in pixels)
  int _screenY;

  /// Display origin (in mm)
  double _origin_x;

  /// Display origin (in mm)
  double _origin_y;

  /// Display origin (in mm)
  double _origin_z;

  /// Display resolution
  double _resolution;

  /// Display axes
  double _xaxis[3], _yaxis[3], _zaxis[3];

  /// Region of interest
  double _x1, _y1, _z1, _x2, _y2, _z2;

  /// Current mouse
  /// Interpolator for target image
  mirtk::InterpolateImageFunction *_targetInterpolator;

  /// Interpolator for source image
  mirtk::InterpolateImageFunction *_sourceInterpolator;

  /// Interpolator for segmentation image
  mirtk::InterpolateImageFunction *_segmentationInterpolator;

  /// Interpolator for selection image
  mirtk::InterpolateImageFunction *_selectionInterpolator;

  /// Flag whether transformation for reslicing of source image should be applied
  bool _sourceTransformApply;

  /// Flag whether transformation for reslicing of source image should be inverted
  bool _sourceTransformInvert;

  /// Display viewing mix in shutter viewing mode
  double _viewMix;

  /// Flag for rview mode
  RViewMode _viewMode;

  /// Flag for configuration mode
  ConfigViewerMode _configMode;

  /// Flag for display labels of segmentation image
  int _DisplaySegmentationLabels;

  /// Flag for display contours of segmentation image
  int _DisplaySegmentationContours;

  /// Flag for segmentation mode
  int _SegmentationMode;

  /// Parameter for paintbrush width
  int _PaintBrushWidth;

  /// Minimum region growing threshold
  int _RegionGrowingThresholdMin;

  /// Maximum region growing threshold
  int _RegionGrowingThresholdMax;

  /// Deformation property
  DeformationProperty _DeformationProperty;

  /// Deformation blending
  double _DeformationBlending;

  /// Flag for display orientation
  DisplayMode _DisplayMode;

  /// Flag for line thickness
  double _LineThickness;

  /// Flag for spedd
  double _Speed;

  /// Flag for display of isolines from target image
  int _DisplayTargetContour;

  /// Flag for display of isolines from source image
  int _DisplaySourceContour;

  /// Flag for snap to grid
  int _SnapToGrid;

  /// Flag for display of cross hair
  int _DisplayCursor;

  /// Mode for cursor display
  CursorMode _CursorMode;

  /// Flag for display of axis labels
  int _DisplayAxisLabels;

  /// Resolution for display of deformation grid
  int _DisplayDeformationGridResolution;

  /// Flag for display of deformation grid
  int _DisplayDeformationGrid;

  /// Flag for display of deformation points
  int _DisplayDeformationPoints;

  /// Flag for display of deformation arrows
  int _DisplayDeformationArrows;

  /// 0: local, 1: global+local
  int _DisplayDeformationTotal;

  /// Flag for display of landmarks
  int _DisplayLandmarks;

  /// IDs of target landmarks to display
  std::set<int> _selectedTargetLandmarks;

  /// IDs of source landmarks to display
  std::set<int> _selectedSourceLandmarks;

  /// Flag for display of ROI
  int _DisplayROI;

  /// Flag for track of tag using gravity window
  int _TrackTAG;

  /// Flag for display of tag grid pattern
  int _ViewTAG;

  /// Flag for flipping X coordinates
  int _FlipX;

  /// Flag for flipping Y coordinates
  int _FlipY;

  /// Flag for flipping Z coordinates
  int _FlipZ;

  /// Current mouse position
  int _mouseX, _mouseY, _mouseZ;

  /// Current intensity at mouse position
  double _mouseTargetIntensity;

  /// Current viewer in which the mouse is
  int _mouseViewer;

  /// Region growing mode
  RegionGrowingMode _regionGrowingMode;

#if MIRTK_IO_WITH_VTK && defined(HAVE_VTK)
  /// Flag for display of object
  int _DisplayObject;

  /// Flag for warping of object
  int _DisplayObjectWarp;

  /// Flag for display of object grid
  int _DisplayObjectGrid;
#endif

public:

  /// Constructor
  RView(int, int);

  /// Destructor
  virtual ~RView();

  /// Render
  void Draw();

  /// Render offscreen
  void DrawOffscreen(char *);

  /// Update registration viewer
  void Update();

  /// Set update of source transformation to on
  void SourceUpdateOn();

  /// Set update of segmentation transformation to on
  void SegmentationUpdateOn();

  /// Resize registration viewer
  void Resize(int, int);

  /// Initialize registration viewer
  virtual void Initialize(bool = true);

  /// Configure registration viewer
  virtual void Configure(RViewConfig []);

  /// Read configuration
  virtual void Read(char *);

  /// Write configuration
  virtual void Write(char *);

  /// Read target image
  virtual void ReadTarget(char *);

  /// Read target image sequence
  virtual void ReadTarget(int, char **);

  /// Read source image
  virtual void ReadSource(char *);

  /// Read source image sequence
  virtual void ReadSource(int, char **);

  /// Read segmentation image
  virtual void ReadSegmentation(char *);

  /// Write target image
  virtual void WriteTarget(char *);

  /// Write source image
  virtual void WriteSource(char *);

  /// Write segmentation image
  virtual void WriteSegmentation(char *);

  /// Read transformation
  virtual void ReadTransformation(char *);

  /// Write transformation
  virtual void WriteTransformation(char *);

  /// Read target landmarks
  virtual void ReadTargetLandmarks(char *);

  /// Read source landmarks
  virtual void ReadSourceLandmarks(char *);

  /// Write target landmarks
  virtual void WriteTargetLandmarks(char *);

  /// Write sourcelandmarks
  virtual void WriteSourceLandmarks(char *);

#if MIRTK_IO_WITH_VTK && defined(HAVE_VTK)
  /// Read object
  virtual void ReadObject(const char *);

  /// Remove object
  virtual void RemoveObject();

  /// Return object
  virtual vtkPointSet *GetObject(int);

  /// Turn display of object movie
  void ObjectMovieOn();

  /// Turn display of object movie
  void ObjectMovieOff();

  /// Return display of object movie
  int GetObjectMovie();
#endif

  /// Get width of registration viewer (in pixels)
  int GetWidth();

  /// Get height of registration viewer (in pixels)
  int GetHeight();

  /// Set display resolution (in mm)
  void   SetResolution(double);

  /// Get display resolution (in mm)
  double GetResolution();

  /// Set display origin (in mm)
  void SetOrigin(double, double, double);

  /// Set display origin (in mm) of target (side-by-side view only)
  void SetTargetOrigin(double, double, double);

  /// Set display origin (in mm) of source (side-by-side view only)
  void SetSourceOrigin(double, double, double);

  /// Set display origin (in display pixels)
  void SetOrigin(int, int);

  /// Get display origin (in mm)
  void GetOrigin(double &, double &, double &);

  /// Set target frame
  void SetTargetFrame(int);

  /// Get target frame
  int GetTargetFrame();

  /// Set source frame
  void SetSourceFrame(int);

  /// Get source frame
  int GetSourceFrame();

  /// Get minimum target intensity
  double GetTargetMin();
  
  /// Get maximum target intensity
  double GetTargetMax();
  
  /// Get minimum source intensity
  double GetSourceMin();
  
  /// Get maximum source intensity
  double GetSourceMax();
  
  /// Get minimum subtraction intensity
  double GetSubtractionMin();
  
  /// Get maximum subtraction intensity
  double GetSubtractionMax();
  
  /// Set ROI to default parameters
  void ResetROI();

  // Update left-upper corner of ROI
  void UpdateROI1(int, int);

  // Update right-lower corner of ROI
  void UpdateROI2(int, int);

  /// Add a point to the contour
  void AddContour(int, int, ContourMode mode);

  /// Undo adding last part contour
  void UndoContour();

  /// Delete current contour
  void ClearContour();

  /// Fill contour
  void FillContour(int, int);

  /// Fills an area defined by contour
  void FillArea(int, int);

  /// Region growing for contour
  void RegionGrowContour(int, int);

  /// Set ROI
  void SetROI(double, double, double, double, double, double);

  /// Get ROI
  void GetROI(double &, double &, double &, double &, double &, double &);

  /// Get an information string
  void GetInfoText(char *, char *, char *, char *, char *);

  /// Current mouse position
  void MousePosition(int, int);

  /// Current mouse wheel
  void MouseWheel(int, int, int);

  /// Set glLine thickness
  void SetLineThickness(double value);

  /// Get glLine thickness
  double GetLineThickness();
  
  /// Set speed
  void SetSpeed(double value);

  /// Get speed
  double GetSpeed();

  /// Get an information string about the transformation level
  void GetTransformationText(std::list<char *> &);

  /// Turn iso-contours extracted from target image on
  void DisplayTargetContoursOn();

  /// Turn iso-contours extracted from target image off
  void DisplayTargetContoursOff();

  /// Return display of target iso-contours
  int GetDisplayTargetContours();

  /// Turn iso-contours extracted from source image on
  void DisplaySourceContoursOn();

  /// Turn iso-contours extracted from source image on
  void DisplaySourceContoursOff();

  /// Return display of source iso-contours
  int GetDisplaySourceContours();

  /// Return display mode
  DisplayMode GetDisplayMode();

  /// Set display mode
  void SetDisplayMode(DisplayMode mode);

  /// Turn caching of displacements on
  void CacheDisplacementsOn();

  /// Turn caching of displacements off
  void CacheDisplacementsOff();

  /// Return displacements caching mode
  int GetCacheDisplacements();

  /// Turn snap to grid on
  void SnapToGridOn();

  /// Turn snap to grid off
  void SnapToGridOff();

  /// Return snap to grid
  int GetSnapToGrid();

  /// Turn axis labels on
  void DisplayAxisLabelsOn();

  /// Turn axis labels off
  void DisplayAxisLabelsOff();

  /// Return display of cross hair
  int GetDisplayCursor();

  /// Return cursor mode
  CursorMode GetCursorMode();

  /// Set cursor mode
  void SetCursorMode(CursorMode mode);

  /// Return minimum display intensity of target image
  double GetDisplayMinTarget();

  /// Sets minimum display intensity of target image
  void SetDisplayMinTarget(double);

  /// Return maximum display intensity of target image
  double GetDisplayMaxTarget();

  /// Sets maximum display intensity of target image
  void SetDisplayMaxTarget(double);

  /// Return minimum display intensity of source image
  double GetDisplayMinSource();

  /// Sets minimum display intensity of source image
  void SetDisplayMinSource(double);

  /// Return maximum display intensity of source image
  double GetDisplayMaxSource();

  /// Sets maximum display intensity of source image
  void SetDisplayMaxSource(double);

  /// Return minimum display intensity of subtraction image
  double GetDisplayMinSubtraction();

  /// Sets minimum display intensity of subtraction image
  void SetDisplayMinSubtraction(double);

  /// Return maximum display intensity of subtraction image
  double GetDisplayMaxSubtraction();

  /// Sets maximum display intensity of subtraction image
  void SetDisplayMaxSubtraction(double);

  /// Return minimum display intensity of deformation 
  double GetDisplayMinDeformation();

  /// Sets minimum display intensity of deformation 
  void SetDisplayMinDeformation(double);

  /// Return maximum display intensity of deformation
  double GetDisplayMaxDeformation();

  /// Sets maximum display intensity of deformation
  void SetDisplayMaxDeformation(double);

  /// Sets fraction of deformation to display
  void SetDisplayDeformationBlending(double);

  /// Return fraction of deformation to display
  double GetDisplayDeformationBlending();

  /// Turn cross hair on
  void DisplayCursorOn();

  /// Turn cross hair off
  void DisplayCursorOff();

  /// Turn display of deformation grid on
  void DisplayDeformationGridOn();

  /// Turn display of deformation grid off
  void DisplayDeformationGridOff();

  /// Return display of deformation grid
  int GetDisplayDeformationGrid();

  /// Return display of deformation grid
  int GetDisplayDeformationGridResolution();

  /// Sets display resolution of deformation grid
  void SetDisplayDeformationGridResolution(int);

  /// Turn display of deformation control points on
  void DisplayDeformationPointsOn();

  /// Turn display of deformation control points off
  void DisplayDeformationPointsOff();

  /// Return display of deformation points
  int GetDisplayDeformationPoints();

  /// Turn display of total deformation on
  void DisplayDeformationTotalOn();

  /// Turn display of total deformation off
  void DisplayDeformationTotalOff();

  /// Return display of total deformation
  int GetDisplayDeformationTotal();

  /// Turn display of ROI
  void DisplayROIOn();

  /// Turn display of ROI
  void DisplayROIOff();

  /// Return display of ROI
  int GetDisplayROI();

  /// Turn display of ROI
  void TrackTAGOn();

  /// Turn display of ROI
  void TrackTAGOff();

  /// Return display of ROI
  int GetTrackTAG();

    /// Turn display of ROI
  void ViewTAGOn();

  /// Turn display of ROI
  void ViewTAGOff();

  /// Return display of ROI
  int GetViewTAG();

  /// Turn display of deformation arrows on
  void DisplayDeformationArrowsOn();

  /// Turn display of deformation arrows off
  void DisplayDeformationArrowsOff();

  /// Return display of deformation arrows
  int GetDisplayDeformationArrows();

  /// Return display of relative deformations
  int GetDisplayDeformationRelative();

  /// Turn display of landmarks on
  void DisplayLandmarksOn();

  /// Turn display of landmarks off
  void DisplayLandmarksOff();

  /// Return display of landmarks
  int GetDisplayLandmarks();

  /// Turn display of specified target landmark on
  void SelectTargetLandmark(int);

  /// Turn display of specified target landmark off
  void DeselectTargetLandmark(int);

  /// Return display of specified target landmark
  int IsTargetLandmarkSelected(int);

  /// Clear selection of target landmarks to display
  void ClearTargetLandmarkSelection();

  /// Turn display of specified source landmark on
  void SelectSourceLandmark(int);

  /// Turn display of specified source landmark off
  void DeselectSourceLandmark(int);

  /// Return display of specified source landmark
  int IsSourceLandmarkSelected(int);

  /// Clear selection of target landmarks to display
  void ClearSourceLandmarkSelection();

#if MIRTK_IO_WITH_VTK && defined(HAVE_VTK)
  /// Turn display of object on
  void DisplayObjectOn();

  /// Turn display of object off
  void DisplayObjectOff();

  /// Return display of object
  int GetDisplayObject();

  /// Turn warping of object on
  void DisplayObjectWarpOn();

  /// Turn warping of object off
  void DisplayObjectWarpOff();

  /// Return warping of object
  int GetDisplayObjectWarp();

  /// Turn grid display of object on
  void DisplayObjectGridOn();

  /// Turn grid display of object off
  void DisplayObjectGridOff();

  /// Return grid display of object
  int GetDisplayObjectGrid();
#endif

  /// Set region growing mode
  void SetRegionGrowingMode(RegionGrowingMode);

  /// Get region growing mode
  RegionGrowingMode GetRegionGrowingMode();

  /// Flip X on
  void FlipXOn();

  /// Flip X off
  void FlipXOff();

  /// Flip Y on
  void FlipYOn();

  /// Flip Y off
  void FlipYOff();

  /// Flip Z on
  void FlipZOn();

  /// Flip Z off
  void FlipZOff();

  /// Set the viewing mix for target and source in shutter viewing mode
  void   SetViewMix(double);

  /// Get the viewing mix for target and source in shutter viewing mode
  double GetViewMix();

  /// Get viewing mode for registration viewer
  RViewMode GetViewMode();

  /// Set viewing mode for registration viewer
  void SetViewMode(RViewMode);

  /// Get configuration mode for viewer
  ConfigViewerMode GetConfigMode();

  /// Set configuration mode for viewer
  void SetConfigMode(ConfigViewerMode );

  /// Set interpolation mode for target image
  void SetTargetInterpolationMode(mirtk::InterpolationMode);

  /// Get interpolation mode for target image
  mirtk::InterpolationMode GetTargetInterpolationMode();

  /// Set interpolation mode for source image
  void SetSourceInterpolationMode(mirtk::InterpolationMode);

  /// Get interpolation model fo target image
  mirtk::InterpolationMode GetSourceInterpolationMode();

  /// Set transformation apply flag for source image
  void SetSourceTransformApply(bool);

  /// Get transformation apply flag for source image
  bool GetSourceTransformApply();

  /// Set transformation invert flag for source image
  void SetSourceTransformInvert(bool);

  /// Get transformation invert flag for source image
  bool GetSourceTransformInvert();

  /// Get a pointer to target image
  mirtk::Image *GetTarget();

  /// Get a pointer to source image
  mirtk::Image *GetSource();

  /// Get a pointer to the lookup table of the target image
  LookupTable *GetTargetLookupTable();

  /// Get a pointer to the lookup table of the source image
  LookupTable *GetSourceLookupTable();

  /// Get a pointer to the lookup table of the subtraction of target and source
  LookupTable *GetSubtractionLookupTable();

  /// Get a pointer to the lookup table of the deformation
  LookupTable *GetDeformationLookupTable();

  /// Get transformation
  mirtk::Transformation *GetTransformation();

  /// Get local transformation
  mirtk::MultiLevelFreeFormTransformation *GetMFFD();

  /// Reset the display origin to origin of target image
  void Reset();

  /// Get a pointer to segment table
  SegmentTable *GetSegmentTable();

  /// Set interpolation mode for segmentation image
  void SetSegmentationInterpolationMode(mirtk::InterpolationMode);

  /// Get interpolation mode for segmentation image
  mirtk::InterpolationMode GetSegmentationInterpolationMode();

  /// Get a pointer to segmentation image
  mirtk::GreyImage *GetSegmentation();

  /// Get a pointer to the lookup table of the segmentation image
  LookupTable *GetSegmentationLookupTable();

  /// Get a pointer to the lookup table of the segmentation image
  VoxelContour *GetVoxelContour();

  /// Turns on the segmentation drawing
  void SegmentationLabelsOn();

  /// Turns off the segmentation drawing
  void SegmentationLabelsOff();

  /// Gets the current display value
  int GetDisplaySegmentationLabels();

  /// Turns on the segmentation drawing
  void SegmentationContoursOn();

  /// Turns off the segmentation drawing
  void SegmentationContoursOff();

  /// Gets the current display value
  int GetDisplaySegmentationContours();

  /// Reset the display origin to origin of segmentation image
  void ResetSegmentation();

  /// Set segmentation mode
  void SegmentationMode(int mode);

  /// Set paint brush width
  void SetPaintBrushWidth(int width);

  /// Set minimum region growing threshold
  void SetRegionGrowingThresholdMinimum(int threshold);

  /// Set maximum region growing threshold
  void SetRegionGrowingThresholdMaximum(int threshold);

  /// Returns segmentation mode
  int GetSegmentationMode();

  /// Returns PaintBrushWidth
  int GetPaintBrushWidth();

  /// Get minimum region growing threshold
  int GetRegionGrowingThresholdMinimum();

  /// Get maximum region growing threshold
  int GetRegionGrowingThresholdMaximum();

  /// Clipping of a drawable to the viewport
  void Clip();

  /// Add target landmark
  void AddTargetLandmark(mirtk::Point &, char *);

  /// Add source landmark
  void AddSourceLandmark(mirtk::Point &, char *);

  /// Delete target landmark
  void DeleteTargetLandmark(int);

  /// Delete source landmark
  void DeleteSourceLandmark(int);

  /// Insert target landmark
  void InsertTargetLandmark(mirtk::Point &, int, char *);

  /// Insert source landmark
  void InsertSourceLandmark(mirtk::Point &, int, char *);

  /// Get target landmark
  void GetTargetLandmark(mirtk::Point &, int, char *);

  /// Get source landmark
  void GetSourceLandmark(mirtk::Point &, int, char *);

  /// Put target landmark
  void PutTargetLandmark(mirtk::Point, int, char *);

  /// Put source landmark
  void PutSourceLandmark(mirtk::Point, int, char *);

  /// Label target landmark
  void LabelTargetLandmark(int, char *);

  /// Label source landmark
  void LabelSourceLandmark(int, char *);

  /// Return number of target landmarks
  int GetNumberOfTargetLandmarks();

  /// Return number of source landmarks
  int GetNumberOfSourceLandmarks();

  /// Approximate transformation using landmark correspondences
  double FitLandmarks();

  /// Callback method for function keys
  void cb_special(int key, int x, int y,
                  int target_delta, int source_delta);

  /// Callback method for function key information
  void cb_special_info();

  /// Callback method for keyboard events
  void cb_keyboard(unsigned char key);

  /// Callback method for keyboard event information
  void cb_keyboard_info();

};

inline void RView::SourceUpdateOn()
{
  _sourceUpdate = true;
}

inline void RView::SegmentationUpdateOn()
{
  _segmentationUpdate = true;
}


inline int RView::GetWidth()
{
  return _screenX;
}

inline int RView::GetHeight()
{
  return _screenY;
}

inline double RView::GetResolution()
{
  return _resolution;
}

inline void RView::SetResolution(double resolution)
{
  _resolution = resolution;
  this->Initialize(false);
}

inline void RView::SetViewMode(RViewMode value)
{
  _viewMode = value;
}

inline RViewMode RView::GetViewMode()
{
  return _viewMode;
}


inline ConfigViewerMode RView::GetConfigMode()
{
  return _configMode;
}

inline void RView::SetConfigMode(ConfigViewerMode configMode)
{
  _configMode = configMode;
}

inline void RView::SetViewMix(double value)
{
  _viewMix = value;
}

inline double RView::GetViewMix()
{
  return _viewMix;
}

inline void RView::SetLineThickness(double value)
{
  _LineThickness = value;
}

inline double RView::GetLineThickness()
{
  return _LineThickness;
}

inline void RView::CacheDisplacementsOn()
{
  _CacheDisplacements = true;
  this->Initialize(true);
}

inline void RView::CacheDisplacementsOff()
{
  _CacheDisplacements = false;
  this->Initialize(true);
}

inline int RView::GetCacheDisplacements()
{
  return _CacheDisplacements;
}

inline void RView::SetSpeed(double value)
{
  _Speed = value;
}

inline double RView::GetSpeed()
{
    if (_Speed >= 0)
        return _Speed;
    else
        return -1/_Speed;
}

inline void RView::DisplayTargetContoursOn()
{
  _DisplayTargetContour = true;
}

inline void RView::DisplayTargetContoursOff()
{
  _DisplayTargetContour = false;
}

inline int RView::GetDisplayTargetContours()
{
  return _DisplayTargetContour;
}

inline void RView::DisplaySourceContoursOn()
{
  _DisplaySourceContour = true;
}

inline void RView::DisplaySourceContoursOff()
{
  _DisplaySourceContour = false;
}

inline int RView::GetDisplaySourceContours()
{
  return _DisplaySourceContour;
}

inline void RView::SnapToGridOn()
{
  _SnapToGrid = true;

  // Round origin to nearest voxel
  _targetImage->WorldToImage(_origin_x, _origin_y, _origin_z);
  _origin_x = round(_origin_x);
  _origin_y = round(_origin_y);
  _origin_z = round(_origin_z);
  _targetImage->ImageToWorld(_origin_x, _origin_y, _origin_z);

  // Update viewer
  for (int k = 0; k < _NoOfViewers; k++) {
    _targetImageOutput[k]->PutOrigin(_origin_x, _origin_y, _origin_z);
    _sourceImageOutput[k]->PutOrigin(_origin_x, _origin_y, _origin_z);
    _segmentationImageOutput[k]->PutOrigin(_origin_x, _origin_y, _origin_z);
    _selectionImageOutput[k]->PutOrigin(_origin_x, _origin_y, _origin_z);
  }

  // Update everything else
  _targetUpdate = true;
  _sourceUpdate = true;
  _segmentationUpdate = true;
  _selectionUpdate = true;

}

inline void RView::SnapToGridOff()
{
  _SnapToGrid = false;
}

inline int RView::GetSnapToGrid()
{
  return _SnapToGrid;
}

inline void RView::DisplayCursorOn()
{
  _DisplayCursor = true;
}

inline void RView::DisplayCursorOff()
{
  _DisplayCursor = false;
}

inline CursorMode RView::GetCursorMode()
{
  return _CursorMode;
}

inline void RView::SetCursorMode(CursorMode mode)
{
  _CursorMode = mode;
}

inline int RView::GetDisplayCursor()
{
  return _DisplayCursor;
}

inline void RView::DisplayAxisLabelsOn()
{
  _DisplayAxisLabels = true;
}

inline void RView::DisplayAxisLabelsOff()
{
  _DisplayAxisLabels = false;
}

inline DisplayMode RView::GetDisplayMode()
{
  return _DisplayMode;
}

inline void RView::SetDisplayMode(DisplayMode mode)
{
  _DisplayMode = mode;
}

inline void RView::DisplayDeformationGridOn()
{
  _DisplayDeformationGrid = true;
}

inline void RView::DisplayDeformationGridOff()
{
  _DisplayDeformationGrid = false;
}

inline int RView::GetDisplayDeformationGrid()
{
  return _DisplayDeformationGrid;
}

inline int RView::GetDisplayDeformationGridResolution()
{
  return _DisplayDeformationGridResolution;
}

inline void RView::SetDisplayDeformationBlending(double a)
{
  if      (a < 0.0) _DeformationBlending = 0.0;
  else if (a > 1.0) _DeformationBlending = 1.0;
  else              _DeformationBlending = a;
  //if (_sourceTransformApply) _sourceUpdate = true;
}

inline double RView::GetDisplayDeformationBlending()
{
  return _DeformationBlending;
}

inline void RView::SetDisplayDeformationGridResolution(int res)
{
  _DisplayDeformationGridResolution = res;
}

inline void RView::DisplayDeformationPointsOn()
{
  _DisplayDeformationPoints = true;
}

inline void RView::DisplayDeformationPointsOff()
{
  _DisplayDeformationPoints = false;
}

inline int RView::GetDisplayDeformationPoints()
{
  return _DisplayDeformationPoints;
}

inline void RView::DisplayDeformationArrowsOn()
{
  _DisplayDeformationArrows = true;
}

inline void RView::DisplayDeformationArrowsOff()
{
  _DisplayDeformationArrows = false;
}

inline int RView::GetDisplayDeformationArrows()
{
  return _DisplayDeformationArrows;
}

inline void RView::DisplayDeformationTotalOn()
{
  _DisplayDeformationTotal = true;
}

inline void RView::DisplayDeformationTotalOff()
{
  _DisplayDeformationTotal = false;
}

inline int RView::GetDisplayDeformationTotal()
{
  return _DisplayDeformationTotal;
}

inline void RView::DisplayLandmarksOn()
{
  _DisplayLandmarks = true;
}

inline void RView::DisplayLandmarksOff()
{
  _DisplayLandmarks = false;
}

inline int RView::GetDisplayLandmarks()
{
  return _DisplayLandmarks;
}

inline void RView::SelectTargetLandmark(int id)
{
  if (0 < id && id <= _targetLandmarks.Size()) {
    _selectedTargetLandmarks.insert(id-1);
  }
}

inline void RView::DeselectTargetLandmark(int id)
{
  _selectedTargetLandmarks.erase(id-1);
}

inline int RView::IsTargetLandmarkSelected(int id)
{
  return _selectedTargetLandmarks.find(id-1) != _selectedTargetLandmarks.end();
}

inline void RView::ClearTargetLandmarkSelection()
{
  _selectedTargetLandmarks.clear();
}

inline void RView::SelectSourceLandmark(int id)
{
  if (0 < id && id <= _targetLandmarks.Size()) {
    _selectedSourceLandmarks.insert(id-1);
  }
}

inline void RView::DeselectSourceLandmark(int id)
{
  _selectedSourceLandmarks.erase(id-1);
}

inline int RView::IsSourceLandmarkSelected(int id)
{
  return _selectedSourceLandmarks.find(id-1) != _selectedSourceLandmarks.end();
}

inline void RView::ClearSourceLandmarkSelection()
{
  _selectedSourceLandmarks.clear();
}

inline int RView::GetDisplaySegmentationLabels()
{
  return _DisplaySegmentationLabels;
}

inline int RView::GetDisplaySegmentationContours()
{
  return _DisplaySegmentationContours;
}

inline int RView::GetSegmentationMode()
{
  return _SegmentationMode;
}

inline int RView::GetPaintBrushWidth()
{
  return _PaintBrushWidth;
}

inline int RView::GetRegionGrowingThresholdMinimum()
{
  return _RegionGrowingThresholdMin;
}

inline int RView::GetRegionGrowingThresholdMaximum()
{
  return _RegionGrowingThresholdMax;
}

inline void RView::DisplayROIOn()
{
  _DisplayROI = true;
}

inline void RView::DisplayROIOff()
{
  _DisplayROI = false;
}

inline int RView::GetDisplayROI()
{
  return _DisplayROI;
}

inline void RView::TrackTAGOn()
{
  _TrackTAG = true;
}

inline void RView::TrackTAGOff()
{
  _TrackTAG = false;
}

inline int RView::GetTrackTAG()
{
  return _TrackTAG;
}

inline void RView::ViewTAGOn()
{
  _ViewTAG = true;
}

inline void RView::ViewTAGOff()
{
  _ViewTAG = false;
}

inline int RView::GetViewTAG()
{
  return _ViewTAG;
}

#if MIRTK_IO_WITH_VTK && defined(HAVE_VTK)

inline vtkPointSet *RView::GetObject(int i)
{
  if ((i < 0) || (i > _NoOfObjects-1)) {
    std::cerr << "RView::GetObject: Invalid object: " << i << std::endl;
    return NULL;
  }
  return _Object[i];
}

inline void RView::ObjectMovieOn()
{
    _ObjectMovie = true;
}

inline void RView::ObjectMovieOff()
{
    _ObjectMovie = false;
}

inline int RView::GetObjectMovie()
{
    return _ObjectMovie;
}

inline void RView::DisplayObjectOn()
{
    _DisplayObject = true;
}

inline void RView::DisplayObjectOff()
{
  _DisplayObject = false;
}

inline int RView::GetDisplayObject()
{
  return _DisplayObject;
}

inline void RView::DisplayObjectWarpOn()
{
  _DisplayObjectWarp = true;
}

inline void RView::DisplayObjectWarpOff()
{
  _DisplayObjectWarp = false;
}

inline int RView::GetDisplayObjectWarp()
{
  return _DisplayObjectWarp;
}

inline void RView::DisplayObjectGridOn()
{
  _DisplayObjectGrid = true;
}

inline void RView::DisplayObjectGridOff()
{
  _DisplayObjectGrid = false;
}

inline int RView::GetDisplayObjectGrid()
{
  return _DisplayObjectGrid;
}

#endif

inline RegionGrowingMode RView::GetRegionGrowingMode()
{
  return _regionGrowingMode;
}

inline void RView::SetRegionGrowingMode(RegionGrowingMode mode)
{
  _regionGrowingMode = mode;
}

inline void RView::FlipXOff()
{
  if (_FlipX == true) {
    _FlipX = false;
    _xaxis[0] *= -1;
    _xaxis[1] *= -1;
    _xaxis[2] *= -1;
    this->Initialize();
  }
}

inline void RView::FlipXOn()
{
  if (_FlipX == false) {
    _FlipX = true;
    _xaxis[0] *= -1;
    _xaxis[1] *= -1;
    _xaxis[2] *= -1;
    this->Initialize();
  }
}

inline void RView::FlipYOff()
{
  if (_FlipY == true) {
    _FlipY = false;
    _yaxis[0] *= -1;
    _yaxis[1] *= -1;
    _yaxis[2] *= -1;
    this->Initialize();
  }
}

inline void RView::FlipYOn()
{
  if (_FlipY == false) {
    _FlipY = true;
    _yaxis[0] *= -1;
    _yaxis[1] *= -1;
    _yaxis[2] *= -1;
    this->Initialize();
  }
}

inline void RView::FlipZOff()
{
  if (_FlipZ == true) {
    _FlipZ = false;
    _zaxis[0] *= -1;
    _zaxis[1] *= -1;
    _zaxis[2] *= -1;
    this->Initialize();
  }
}

inline void RView::FlipZOn()
{
  if (_FlipZ == false) {
    _FlipZ = true;
    _zaxis[0] *= -1;
    _zaxis[1] *= -1;
    _zaxis[2] *= -1;
    this->Initialize();
  }
}

inline void RView::SetOrigin(double x, double y, double z)
{
  int i;

  _origin_x = x;
  _origin_y = y;
  _origin_z = z;
  for (i = 0; i < _NoOfViewers; i++) {
    _targetImageOutput[i]->PutOrigin(_origin_x, _origin_y, _origin_z);
    _sourceImageOutput[i]->PutOrigin(_origin_x, _origin_y, _origin_z);
    _segmentationImageOutput[i]->PutOrigin(_origin_x, _origin_y, _origin_z);
    _selectionImageOutput[i]->PutOrigin(_origin_x, _origin_y, _origin_z);
  }
  _targetUpdate       = true;
  _sourceUpdate       = true;
  _segmentationUpdate = true;
  _selectionUpdate    = true;
}

inline void RView::SetTargetOrigin(double x, double y, double z)
{
  _origin_x = x;
  _origin_y = y;
  _origin_z = z;
  for (int i = 0; i < _NoOfViewers; ++i) {
    if (!_isSourceViewer[i]) {
      _targetImageOutput      [i]->PutOrigin(_origin_x, _origin_y, _origin_z);
      _sourceImageOutput      [i]->PutOrigin(_origin_x, _origin_y, _origin_z);
      _segmentationImageOutput[i]->PutOrigin(_origin_x, _origin_y, _origin_z);
      _selectionImageOutput   [i]->PutOrigin(_origin_x, _origin_y, _origin_z);
    }
  }
  _targetUpdate       = true;
  _sourceUpdate       = true;
  _segmentationUpdate = true;
  _selectionUpdate    = true;
}

inline void RView::SetSourceOrigin(double x, double y, double z)
{
  for (int i = 0; i < _NoOfViewers; ++i) {
    if (_isSourceViewer[i]) {
      _targetImageOutput      [i]->PutOrigin(x, y, z);
      _sourceImageOutput      [i]->PutOrigin(x, y, z);
      _segmentationImageOutput[i]->PutOrigin(x, y, z);
      _selectionImageOutput   [i]->PutOrigin(x, y, z);
    }
  }
  _targetUpdate       = true;
  _sourceUpdate       = true;
  _segmentationUpdate = true;
  _selectionUpdate    = true;
}

inline void RView::GetOrigin(double &x, double &y, double &z)
{
  x = _origin_x;
  y = _origin_y;
  z = _origin_z;
}

inline void RView::SetROI(double x1, double y1, double z1, double x2, double y2, double z2)
{
  _x1 = x1;
  _y1 = y1;
  _z1 = z1;
  _x2 = x2;
  _y2 = y2;
  _z2 = z2;

}

inline void RView::GetROI(double &x1, double &y1, double &z1, double &x2, double &y2, double &z2)
{
  x1 = _x1;
  y1 = _y1;
  z1 = _z1;
  x2 = _x2;
  y2 = _y2;
  z2 = _z2;
}

inline mirtk::Image *RView::GetTarget()
{
  return _targetImage;
}

inline mirtk::Image *RView::GetSource()
{
  return _sourceImage;
}

inline mirtk::GreyImage *RView::GetSegmentation()
{
  return _segmentationImage;
}

inline VoxelContour *RView::GetVoxelContour()
{
  return &_voxelContour;
}

inline SegmentTable *RView::GetSegmentTable()
{
  return _segmentTable;
}

inline LookupTable *RView::GetTargetLookupTable()
{
  return _targetLookupTable;
}

inline LookupTable *RView::GetSourceLookupTable()
{
  return _sourceLookupTable;
}

inline LookupTable *RView::GetSubtractionLookupTable()
{
  return _subtractionLookupTable;
}

inline mirtk::Transformation *RView::GetTransformation()
{
  return _sourceTransform;
}

inline mirtk::MultiLevelFreeFormTransformation *RView::GetMFFD()
{
  mirtk::MultiLevelFreeFormTransformation *transform = dynamic_cast<mirtk::MultiLevelFreeFormTransformation *>(_sourceTransform);
  return transform;
}

inline void RView::Clip()
{
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
  glViewport(0, 0, (GLsizei) _screenX, (GLsizei) _screenY);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(0.0, (GLdouble) _screenX, 0.0, (GLdouble) _screenY);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
}

inline void RView::AddTargetLandmark(mirtk::Point &point, char *)
{
  // Add landmark as point, ignoring label for now
      _targetLandmarks.Add(point);
}

inline void RView::AddSourceLandmark(mirtk::Point &point, char *)
{
  // Add landmark as point, ignoring label for now	
	  _sourceLandmarks.Add(point);
}

inline void RView::DeleteTargetLandmark(int id)
{
  // Delete landmark from list
  if ((id > 0) && (id <= _targetLandmarks.Size())) {
    DeselectTargetLandmark(id);
    mirtk::Point p = _targetLandmarks(id-1);
    _targetLandmarks.Del(p);
  }
}

inline void RView::DeleteSourceLandmark(int id)
{
  // Delete landmark from list
  if ((id > 0) && (id <= _sourceLandmarks.Size())) {
    DeselectSourceLandmark(id);
    mirtk::Point p = _sourceLandmarks(id-1);
    _sourceLandmarks.Del(p);
  }
}

inline void RView::InsertTargetLandmark(mirtk::Point &point, int id, char *)
{
  // Insert landmark, ignoring label for now
  if (_targetLandmarks.Size() == 0) {
    _targetLandmarks.Add(point);
  } else if ((id > 0) && (id <= _targetLandmarks.Size())) {
    // Would be nice to call an mirtk::PointSet::Insert method...
    int i;
    mirtk::PointSet pset(_targetLandmarks);
    _targetLandmarks.Clear();
    for (i = 1; i < id; i++) {
      _targetLandmarks.Add(pset(i-1));
    }
    _targetLandmarks.Add(point);
    for (i = id; i <= pset.Size(); i++) {
      _targetLandmarks.Add(pset(i-1));
    }
  } else {
    std::cerr << "RView::InsertTargetLandmark : invalid position " << id << std::endl;
  }
}

inline void RView::InsertSourceLandmark(mirtk::Point &point, int id, char *)
{
  // Insert landmark, ignoring label for now
  if (_sourceLandmarks.Size() == 0) {
    _sourceLandmarks.Add(point);
  } else if ((id > 0) && (id <= _sourceLandmarks.Size())) {
    // Would be nice to call an mirtk::PointSet::Insert method...
    int i;
    mirtk::PointSet pset(_sourceLandmarks);
    _sourceLandmarks.Clear();
    for (i = 1; i < id; i++) {
      _sourceLandmarks.Add(pset(i-1));
    }
    _sourceLandmarks.Add(point);
    for (i = id; i <= pset.Size(); i++) {
      _sourceLandmarks.Add(pset(i-1));
    }
  } else {
    std::cerr << "RView::InsertSourceLandmark : invalid position " << id << std::endl;
  }
}

inline void RView::LabelTargetLandmark(int, char *)
{
  // So far, labelling of mirtk::PointSet is not possible
}

inline void RView::LabelSourceLandmark(int, char *)
{
  // So far, labelling of mirtk::PointSet is not possible
}

inline void RView::GetTargetLandmark(mirtk::Point &point, int id, char *)
{
  // Get landmark from list, ignoring label for now
  if ((id > 0) && (id <= _targetLandmarks.Size())) {
    point = _targetLandmarks(id-1);
  }
}

inline void RView::GetSourceLandmark(mirtk::Point &point, int id, char *)
{
  // Get landmark from list, ignoring label for now
  if ((id > 0) && (id <= _sourceLandmarks.Size())) {
    point = _sourceLandmarks(id-1);
  }
}

inline void RView::PutTargetLandmark(mirtk::Point point, int id, char *)
{
  // Put landmark in list, ignoring label for now
  if ((id > 0) && (id <= _targetLandmarks.Size())) {
    _targetLandmarks(id-1) = point;
  }
}

inline void RView::PutSourceLandmark(mirtk::Point point, int id, char *)
{
  // Put landmark in list, ignoring label for now
  if ((id > 0) && (id <= _sourceLandmarks.Size())) {
    _sourceLandmarks(id-1) = point;
  }
}

inline int RView::GetNumberOfTargetLandmarks()
{
  return _targetLandmarks.Size();
}

inline int RView::GetNumberOfSourceLandmarks()
{
  return _sourceLandmarks.Size();
}

inline void RView::SegmentationLabelsOn()
{
  _DisplaySegmentationLabels = true;
}

inline void RView::SegmentationLabelsOff()
{
  _DisplaySegmentationLabels = false;
}

inline void RView::SegmentationContoursOn()
{
  _DisplaySegmentationContours = true;
}

inline void RView::SegmentationContoursOff()
{
  _DisplaySegmentationContours = false;
}

inline double RView::GetTargetMin()
{
  return _targetMin;
}

inline double RView::GetTargetMax()
{
  return _targetMax;
}

inline double RView::GetSourceMin()
{
  return _sourceMin;
}

inline double RView::GetSourceMax()
{
  return _sourceMax;
}

inline double RView::GetSubtractionMin()
{
  return _subtractionMin;
}

inline double RView::GetSubtractionMax()
{
  return _subtractionMax;
}

inline double RView::GetDisplayMinTarget()
{
	return _targetDisplayMin;
}

inline double RView::GetDisplayMaxTarget()
{
	return _targetDisplayMax;
}

inline void RView::SetDisplayMinTarget(double value)
{
	_targetDisplayMin = value;
	_targetLookupTable->SetMinDisplayIntensity(round((value - _targetMin) * 10000.0 / (_targetMax - _targetMin)));
}

inline void RView::SetDisplayMaxTarget(double value)
{
	_targetDisplayMax = value;
	_targetLookupTable->SetMaxDisplayIntensity(round((value - _targetMin) * 10000.0 / (_targetMax - _targetMin)));
}

inline double RView::GetDisplayMinSource()
{
	return _sourceDisplayMin;
}

inline double RView::GetDisplayMaxSource()
{
	return _sourceDisplayMax;
}

inline void RView::SetDisplayMinSource(double value)
{
	_sourceDisplayMin = value;
	_sourceLookupTable->SetMinDisplayIntensity(round((value - _sourceMin) * 10000.0 / (_sourceMax - _sourceMin)));
}

inline void RView::SetDisplayMaxSource(double value)
{
	_sourceDisplayMax = value;
	_sourceLookupTable->SetMaxDisplayIntensity(round((value - _sourceMin) * 10000.0 / (_sourceMax - _sourceMin)));
}

inline double RView::GetDisplayMinSubtraction()
{
	return _subtractionDisplayMin;
}

inline double RView::GetDisplayMaxSubtraction()
{
	return _subtractionDisplayMax;
}

inline void RView::SetDisplayMinSubtraction(double value)
{
	_subtractionDisplayMin = value;
	_subtractionLookupTable->SetMinDisplayIntensity(round(value * 20000.0 / (_subtractionMax - _subtractionMin)));
}

inline void RView::SetDisplayMaxSubtraction(double value)
{
	_subtractionDisplayMax = value;
	_subtractionLookupTable->SetMaxDisplayIntensity(round(value * 20000.0 / (_subtractionMax - _subtractionMin)));
}

#endif
