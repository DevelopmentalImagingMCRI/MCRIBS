/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2016 Imperial College London
 * Copyright 2013-2016 Andreas Schuh
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

#ifndef MIRTK_DeformableSurfaceModel_H
#define MIRTK_DeformableSurfaceModel_H

#include "mirtk/ObjectiveFunction.h"

#include "mirtk/Array.h"

#include "mirtk/RegisteredImage.h"
#include "mirtk/RegisteredPointSet.h"

#include "mirtk/EnergyTerm.h"
#include "mirtk/ExternalForce.h"
#include "mirtk/InternalForce.h"
#include "mirtk/TransformationConstraint.h"
#include "mirtk/MeshSmoothing.h"

#include "vtkSmartPointer.h"
#include "vtkPointSet.h"
#include "vtkPolyData.h"
#include "vtkAbstractCellLocator.h"


namespace mirtk {


/**
 * Energy function describing a deformable surface model
 */
class DeformableSurfaceModel : public ObjectiveFunction
{
  mirtkObjectMacro(DeformableSurfaceModel);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Input point set / surface mesh
  mirtkPublicAttributeMacro(vtkSmartPointer<vtkPointSet>, Input);

  /// Intensity image
  mirtkPublicAggregateMacro(RegisteredImage, Image);

  /// Implicit surface distance image
  mirtkPublicAggregateMacro(RegisteredImage, ImplicitSurface);

  /// Transformation to deform the point set with (optional)
  mirtkPublicAggregateMacro(class Transformation, Transformation);

  /// Deformed point set / surface mesh
  mirtkReadOnlyAttributeMacro(RegisteredPointSet, PointSet);

  /// Number of energy terms
  mirtkReadOnlyAttributeMacro(int, NumberOfTerms);

  /// Minimum/default radius of node neighborhood
  mirtkPublicAttributeMacro(int, NeighborhoodRadius);

  /// Number of gradient averaging iterations
  mirtkPublicAttributeMacro(int, GradientAveraging);

  /// Weighting function used for averaging of gradient vectors
  mirtkPublicAttributeMacro(MeshSmoothing::WeightFunction, GradientWeighting);

  /// Whether to only average gradient vectors pointing in the same direction
  /// as the unsmoothed gradient at the central node (i.e., positive dot product)
  mirtkPublicAttributeMacro(bool, AverageSignedGradients);

  /// Whether to only average the magnitude of the gradient vectors
  mirtkPublicAttributeMacro(bool, AverageGradientMagnitude);

  /// Minimum (average) output mesh edge length
  mirtkPublicAttributeMacro(double, MinEdgeLength);

  /// Maximum (average) output mesh edge length
  mirtkPublicAttributeMacro(double, MaxEdgeLength);

  /// Minimum edge end point normal angle of feature edges
  mirtkPublicAttributeMacro(double, MinFeatureAngle);

  /// Maximum edge end point normal angle of feature edges
  mirtkPublicAttributeMacro(double, MaxFeatureAngle);

  /// Whether inversion of triangles is allowed during local remeshing
  mirtkPublicAttributeMacro(bool, AllowTriangleInversion);

  /// Remesh deformed surface every n-th iteration
  mirtkPublicAttributeMacro(int, RemeshInterval);

  /// Number of iterations since last performed remeshing
  mirtkAttributeMacro(int, RemeshCounter);

  /// Remesh surface using an adaptive edge length interval based on local curvature
  mirtkPublicAttributeMacro(bool, RemeshAdaptively);

  /// Low-pass filter surface mesh every n-th iteration
  mirtkPublicAttributeMacro(int, LowPassInterval);

  /// Number of low-pass filtering iterations
  mirtkPublicAttributeMacro(int, LowPassIterations);

  /// Low-pass filter band
  mirtkPublicAttributeMacro(double, LowPassBand);

  /// Maximum distance of deformed surface points from input surface
  mirtkPublicAttributeMacro(double, MaxInputDistance);
  vtkSmartPointer<vtkAbstractCellLocator> _InputCellLocator;

  /// Enforce non-self-intersection of deformed surface mesh
  mirtkPublicAttributeMacro(bool, HardNonSelfIntersection);

  /// Minimum required distance between non-adjacent triangular faces
  ///
  /// This distance threshold applies when the center point of the other
  /// face is in front of the current triangle whose collisions with other
  /// triangles is being determined (cf. SurfaceCollisions::FrontfaceCollision).
  mirtkPublicAttributeMacro(double, MinFrontfaceDistance);

  /// Minimum required distance between non-adjacent triangular faces
  ///
  /// This distance threshold applies when the center point of the other
  /// face is at the backside of the current triangle whose collisions with other
  /// triangles is being determined (cf. SurfaceCollisions::FrontfaceCollision).
  mirtkPublicAttributeMacro(double, MinBackfaceDistance);

  /// Maximum angle between face normal and center to closest point vector
  /// required for collision to be detected
  mirtkPublicAttributeMacro(double, MaxCollisionAngle);

  /// Use fast approximate surface triangle collision test
  mirtkPublicAttributeMacro(bool, FastCollisionTest);

  /// Disallow passive nodes from moving
  mirtkPublicAttributeMacro(bool, FixPassivePoints);

  /// Allow nodes to move in outwards normal direction
  mirtkPublicAttributeMacro(bool, AllowExpansion);

  /// Allow nodes to move in inwards normal direction
  mirtkPublicAttributeMacro(bool, AllowContraction);

  /// Whether the deformable model is a surface mesh
  mirtkReadOnlyAttributeMacro(bool, IsSurfaceMesh);

  /// Whether to only minimize the energy of external forces
  mirtkPublicAttributeMacro(bool, MinimizeExtrinsicEnergy);

protected:

  /// Number of iterations since last low-pass filtering
  mutable int _LowPassCounter;

  /// Energy terms corresponding to external forces
  Array<class ExternalForce *> _ExternalForce;
  Array<bool>                  _ExternalForceOwner;

  /// Energy terms corresponding to internal forces
  Array<class InternalForce *> _InternalForce;
  Array<bool>                  _InternalForceOwner;

  /// Energy terms which regularize the parametric transformation
  Array<TransformationConstraint *> _Constraint;
  Array<bool>                       _ConstraintOwner;

  /// Input surface meshes which the deformed surface mesh may not intersect
  Array<vtkSmartPointer<vtkPolyData> > _BoundaryConstraint;

public:

  /// Output surface mesh
  vtkSmartPointer<vtkPointSet> Output() const;

  // ---------------------------------------------------------------------------
  // Construction/Destruction
private:

  /// Copy constructor
  /// \note Intentionally not implemented!
  DeformableSurfaceModel(const DeformableSurfaceModel &);

  /// Assignment operator
  /// \note Intentionally not implemented!
  DeformableSurfaceModel &operator =(const DeformableSurfaceModel &);

public:

  /// Constructor
  DeformableSurfaceModel();

  /// Destructor
  virtual ~DeformableSurfaceModel();

  // ---------------------------------------------------------------------------
  // Energy terms

  /// Initialize energy terms once input and parameters have been set
  virtual void Initialize();

  /// Delete previously added energy terms
  void Clear();

  /// Whether energy function has no terms
  bool Empty() const;

  /// Number of energy terms
  int NumberOfForces() const;

  /// Number of internal force terms
  int NumberOfInternalForces() const;

  /// Number of external force terms
  int NumberOfExternalForces() const;

  /// Add external force term and take over ownership of the object
  void Add(class ExternalForce *, bool = true);

  /// Remove external force term and revoke ownership of the object
  void Sub(class ExternalForce *);

  /// Add internal force term and take over ownership of the object
  void Add(class InternalForce *, bool = true);

  /// Remove internal force term and revoke ownership of the object
  void Sub(class InternalForce *);

  /// Add transformation regularization term and take over ownership of the object
  void Add(TransformationConstraint *, bool = true);

  /// Remove transformation regularization term and revoke ownership of the object
  void Sub(TransformationConstraint *);

  /// Get the n-th energy term
  EnergyTerm *Term(int);

  /// Get the n-th energy term
  const EnergyTerm *Term(int) const;

  /// Get the n-th external force term
  class ExternalForce *ExternalForce(int);

  /// Get the n-th external force term
  const class ExternalForce *ExternalForce(int) const;

  /// Get the n-th internal force term
  class InternalForce *InternalForce(int);

  /// Get the n-th internal force term
  const class InternalForce *InternalForce(int) const;

  /// Determine whether a given force term is an external force
  static bool IsExternalForce(const EnergyTerm *);

  /// Determine whether a given force term is an external implicit surface force
  static bool IsImplicitSurfaceForce(const EnergyTerm *);

  /// Determine whether a given force term is an internal force
  static bool IsInternalForce(const EnergyTerm *);

  /// Determine whether n-th energy term is an external force
  bool IsExternalForce(int) const;

  /// Determine whether n-th energy term is an external implicit surface force
  bool IsImplicitSurfaceForce(int) const;

  /// Determine whether n-th energy term is an internal force
  bool IsInternalForce(int) const;

  // ---------------------------------------------------------------------------
  // Settings

  // Import other overloads
  using ObjectiveFunction::Parameter;

  /// Set parameter value from string
  virtual bool Set(const char *, const char *);

  /// Get parameter as key/value as string map
  virtual ParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Function parameters

  /// Get number of deformable surface parameters
  virtual int NumberOfDOFs() const;

  /// Get number of deformable surface points
  int NumberOfPoints() const;

  /// Set deformable surface parameters
  ///
  /// This function can be used to restore the deformable surface parameters
  /// after a failed update which did not result in the desired improvement.
  ///
  /// \param[in] x Value of deformable surface parameters (DoFs).
  virtual void Put(const double *x);

  /// Get deformable surface parameters
  ///
  /// This function can be used to store a backup of the current deformable
  /// surface parameters before an update such that these can be restored using
  /// the Put member function if the update did not result in the desired change
  /// of the overall energy.
  ///
  /// \param[in] x Current values of deformable surface parameters (DoFs).
  virtual void Get(double *x) const;

  /// Get function parameter value
  ///
  /// \returns Value of specified function parameter (DoF).
  virtual double Get(int) const;

  /// Add change (i.e., scaled gradient) to each deformable surface parameter
  ///
  /// This function updates each parameter of the deformable surface model
  /// given a vector of desired changes, i.e., the computed gradient of the
  /// energy function.
  ///
  /// \param[in,out] dx Change of each function parameter (DoF) as computed by the
  ///                   Gradient member function and scaled by a chosen step length.
  ///                   The change of a parameter may be modified by this function
  ///                   in order to satisfy the hard constraints (if any).
  ///
  /// \returns Maximum change of transformation parameter.
  virtual double Step(double *dx);

  /// Update internal state after change of parameters
  virtual void Update(bool = true);

  /// Update energy function after convergence
  ///
  /// For example, fiducial registration error (FRE) terms may update the
  /// point correspondences before another gradient-based optimization of
  /// the new FRE term.
  ///
  /// \returns Whether the energy function has changed.
  virtual bool Upgrade();

  /// Perform local adaptive remeshing
  ///
  /// \returns Whether the surface model has been remeshed. The Update function
  ///          must be called after such remeshing operation before evaluating
  ///          the energy and/or gradient of the deformable surface model.
  virtual bool Remesh();

  // ---------------------------------------------------------------------------
  // Evaluation

  /// Query/evaluate initial value of energy function
  double InitialValue();

  /// Get initial value of n-th energy term
  double InitialValue(int);

  /// Evaluate energy function
  double Value();

  /// Get value of n-th energy term computed upon last evaluation
  double Value(int);

  /// Evaluate gradient of energy function
  ///
  /// This gradient corresponds to the weighted sum of external and internal
  /// forces of the deformable surface model.
  ///
  /// \param[in]  dx      Gradient of energy function.
  /// \param[in]  step    Step length for finite differences.
  /// \param[out] sgn_chg Whether function parameter value is allowed to
  ///                     change sign when stepping along the computed gradient.
  void Gradient(double *dx, double step = .0, bool *sgn_chg = NULL);

  /// Compute norm of gradient of energy function
  ///
  /// This norm can, for example, be the maximum absolute parameter change,
  /// the maximum control point displacement if a FFD transformation is used
  /// to deform the initial surface mesh, or the maximum vertex displacement.
  double GradientNorm(const double *) const;

  /// Adjust step length range
  ///
  /// \param[in]    dx  Gradient of objective function.
  /// \param[inout] min Minimum step length.
  /// \param[inout] max Maximum step length.
  void GradientStep(const double *dx, double &min, double &max) const;

  /// Evaluate energy function
  ///
  /// This function first updates the internal state of the function object
  /// if required due to a previous change of the function parameters and then
  /// evaluates the current energy function value.
  ///
  /// \param[in]  step    Step length for finite differences.
  /// \param[out] dx      Gradient of objective function.
  ///                     If \c NULL, only the function value is computed.
  /// \param[out] sgn_chg Whether function parameter value is allowed to
  ///                     change sign when stepping along the computed gradient.
  ///                     Ignord if \p dx is \c NULL.
  virtual double Evaluate(double *dx = NULL, double step = .0, bool *sgn_chg = NULL);

protected:

  /// Smooth gradient such that neighboring points move coherently
  virtual void SmoothGradient(double *dx) const;

  /// Adjust node displacements to avoid self-intersections and collisions
  ///
  /// \param[in,out] dx   (Scaled) gradient of objective function.
  /// \param[in]     nsi  Enforce non-self-intersection.
  /// \param[in]     mind Minimum front-facing distance.
  /// \param[in]     minw Minimum back-facing distance.
  void ResolveSurfaceCollisions(double *dx, bool nsi, double mind, double minw) const;

  /// Enforce hard constraints on surface model deformation
  ///
  /// This function clamps a nodes' displacement vector (velocity times \f$\delta t\f$),
  /// if otherwise the hard constraints of the deformable surface model would
  /// be violated. Common hard constraints are non-self-intersection and a
  /// maximum total node displacement. If the surface model is deformed by
  /// a parametric transformation, this function does nothing as hard constraints
  /// can only be enforced during the optimization when the parameters of the
  /// deformable surface model are the positions of the individual surface nodes.
  ///
  /// \param[inout] dx (Scaled) gradient of objective function.
  virtual void EnforceHardConstraints(double *dx) const;

  // ---------------------------------------------------------------------------
  // Debugging
public:

  /// Get unweighted and unnormalized value of n-th energy term
  /// \remarks Use for progress reporting only.
  double RawValue(int);

  /// Write input of data force terms
  virtual void WriteDataSets(const char *, const char *, bool = true) const;

  /// Write gradient of force terms
  virtual void WriteGradient(const char *, const char *) const;

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
inline vtkSmartPointer<vtkPointSet> DeformableSurfaceModel::Output() const
{
  return _PointSet.PointSet();
}

// -----------------------------------------------------------------------------
inline int DeformableSurfaceModel::NumberOfPoints() const
{
  return _PointSet.NumberOfPoints();
}


} // namespace mirtk

#endif // MIRTK_DeformableSurfaceModel_H
