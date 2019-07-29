/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2016 Imperial College London
 * Copyright 2016 Andreas Schuh
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

#include "mirtk/MeanCurvatureConstraint.h"

#include "mirtk/Math.h"
#include "mirtk/Matrix3x3.h"
#include "mirtk/Triangle.h"

#include "mirtk/MeshSmoothing.h"
#include "mirtk/SurfaceCurvature.h"

#include "mirtk/VtkMath.h"
#include "vtkPointData.h"
#include "vtkDataArray.h"
#include "vtkIdList.h"


#define USE_CURVATURE_WEIGHTED_SPRING_FORCE 1


namespace mirtk {


// Register energy term with object factory during static initialization
mirtkAutoRegisterEnergyTermMacro(MeanCurvatureConstraint);


// =============================================================================
// Auxiliaries
// =============================================================================

namespace MeanCurvatureConstraintUtils {


// -----------------------------------------------------------------------------
/// Get indices of other vertices of triangles sharing an edge
///
/// \param[in]  surface Triangulated surface mesh.
/// \param[in]  i       Index of edge start point.
/// \param[in]  k       Index of edge end point.
/// \param[out] j       Index of other point on the "left"  of the edge.
/// \param[out] l       Index of other point on the "right" of the edge.
///
/// \returns Number of adjacent triangles.
int GetEdgeNeighborPoints(vtkPolyData *surface, int i, int k, int &j, int &l)
{
  j = l = -1;
  unsigned short ncells1, ncells2, ncells = 0;
  vtkIdType      *cells1, *cells2, npts, *pts;
  surface->GetPointCells(static_cast<vtkIdType>(i), ncells1, cells1);
  surface->GetPointCells(static_cast<vtkIdType>(k), ncells2, cells2);
  for (unsigned short idx1 = 0; idx1 < ncells1; ++idx1)
  for (unsigned short idx2 = 0; idx2 < ncells2; ++idx2) {
    if (cells1[idx1] == cells2[idx2]) {
      ++ncells;
      if (ncells < 3) {
        surface->GetCellPoints(cells1[idx1], npts, pts);
        if (npts == 3) {
          if (pts[0] == i) {
            if (pts[1] == k) l = pts[2];
            else             j = pts[1];
          } else if (pts[1] == i) {
            if (pts[2] == k) l = pts[0];
            else             j = pts[2];
          } else {
            if (pts[0] == k) l = pts[1];
            else             j = pts[0];
          }
        }
      }
    }
  }
  if (ncells == 1) {
    if (j == -1 && l == -1) ncells = 0;
  } else if (ncells == 2) {
    if (j == -1) --ncells;
    if (l == -1) --ncells;
  }
  return ncells;
}

// -----------------------------------------------------------------------------
/// Calculate mean curvature (cf. vtkCurvatures::GetMeanCurvature)
struct CalculateMeanCurvature
{
  vtkPolyData     *_Surface;
  const EdgeTable *_EdgeTable;
  vtkDataArray    *_Curvature;
  int             *_Count;

//  void operator ()(const blocked_range<int> &ptIds) const
//  {
//    double     v_i[3], v_j[3], v_k[3], v_l[3];
//    double     l_ik, A_ijk, A_ikl, A_sum, cs, sn, angle, H;
//    Vector3    p_i, e_ij, e_ik, e_il, n_ijk, n_ikl, n_cross;
//    int        numAdjPts, adjPtIdx, i, j, k, l;
//    const int *adjPtIds;
//
//    for (i = ptIds.begin(); i != ptIds.end(); ++i) {
//
//      _Surface->GetPoint(i, v_i);
//      _EdgeTable->GetAdjacentPoints(i, numAdjPts, adjPtIds);
//      for (adjPtIdx = 0; adjPtIdx < numAdjPts; ++adjPtIdx) {
//        k = adjPtIds[adjPtIdx];
//        if (i < k && GetEdgeNeighborPoints(_Surface, i, k, j, l) == 2 && j != -1 && l != -1) {
//
//          // Get points of adjacent triangle vertices
//          _Surface->GetPoint(j, v_j);
//          _Surface->GetPoint(k, v_k);
//          _Surface->GetPoint(l, v_l);
//
//          // Compute required vector quantities from vertex points
//          p_i   = Vector3(v_i);
//          e_ij  = Vector3(v_j), e_ij -= p_i;
//          e_ik  = Vector3(v_k), e_ik -= p_i;
//          e_il  = Vector3(v_l), e_il -= p_i;
//          n_ijk = e_ij.Cross(e_ik);
//          n_ikl = e_ik.Cross(e_il);
//
//          // Compute areas of adjacent triangles
//          // Note: The factor 1/2 is cancelled by the factor 2 of atan2
//          A_ijk = n_ijk.Normalize();
//          A_ikl = n_ikl.Normalize();
//          A_sum = A_ijk + A_ikl;
//
//          // Compute cosine and sine of angle made up by the face normals
//          cs = n_ijk.Dot(n_ikl);
//          sn = n_ikl.Cross(n_ijk).Dot(e_ik);
//          if (sn != 0. || cs != 0.) {
//
//            // Compute length of edge
//            l_ik = e_ik.Length();
//
//            // Compute signed angle made up by normals in [-pi, pi]
//            sn   /= l_ik;
//            angle = atan2(sn, cs);
//
//            // Add mean curvature contribution of this edge (excl. factor 3)
//            H = l_ik * angle;
//            if (A_sum > 1e-9) H /= A_sum;
//            _Curvature->SetComponent(i, 0, _Curvature->GetComponent(i, 0) + H);
//            _Curvature->SetComponent(k, 0, _Curvature->GetComponent(k, 0) + H);
//          }
//        }
//      }
//    }
//  }

  void operator ()(const blocked_range<vtkIdType> &cellIds) const
  {
    vtkIdType i, j, k, l;
    vtkIdType numCellPts, *cellPts, cellPtIdx;
    vtkIdType numNborPts, *nborPts, nborPtIdx;
    double    v_i[3], v_j[3], v_k[3], v_l[3];
    double    A_ijk, A_ikl, A_sum, length, cs, sn, angle, H;
    Vector3   p_i, e_ij, e_ik, e_il, n_ijk, n_ikl, n_cross;

    vtkSmartPointer<vtkIdList> nborCellIds;
    nborCellIds = vtkSmartPointer<vtkIdList>::New();

    for (vtkIdType cellId = cellIds.begin(); cellId != cellIds.end(); ++cellId) {
      _Surface->GetCellPoints(cellId, numCellPts, cellPts);
      for (cellPtIdx = 0; cellPtIdx < numCellPts; ++cellPtIdx) {

        // Get start and end point indices of current cell edge
        i = cellPts[cellPtIdx];
        k = cellPts[(cellPtIdx + 1) % numCellPts];

        // Get neighboring face sharing this edge
        _Surface->GetCellEdgeNeighbors(cellId, i, k, nborCellIds);
        if (nborCellIds->GetNumberOfIds() == 1 && nborCellIds->GetId(0) > cellId) {

          // Get indices of two other points of this and the neighboring cell,
          // while preserving the order, i.e., orientation for both cells
          l = cellPts[(cellPtIdx + 2) % numCellPts];
          _Surface->GetCellPoints(nborCellIds->GetId(0), numNborPts, nborPts);
          if (numNborPts > 2) {
            nborPtIdx = 0;
            while (nborPts[nborPtIdx] != i) ++nborPtIdx;
            j = nborPts[(nborPtIdx + 1) % numNborPts];
            if (j == k) j = nborPts[(nborPtIdx + numNborPts - 1) % numNborPts];

            // Get vertex points
            _Surface->GetPoint(i, v_i);
            _Surface->GetPoint(j, v_j);
            _Surface->GetPoint(k, v_k);
            _Surface->GetPoint(l, v_l);

            // Compute required vector quantities from vertex points
            p_i   = Vector3(v_i);
            e_ij  = Vector3(v_j), e_ij -= p_i;
            e_ik  = Vector3(v_k), e_ik -= p_i;
            e_il  = Vector3(v_l), e_il -= p_i;
            n_ijk = e_ij.Cross(e_ik);
            n_ikl = e_ik.Cross(e_il);

            // Compute areas of adjacent triangles
            // Note: The area factor 1/2 is cancelled by the factor 2 of atan2
            A_ijk = n_ijk.Normalize();
            A_ikl = n_ikl.Normalize();
            A_sum = A_ijk + A_ikl;

            // Compute cosine and sine of angle made up by the face normals
            cs = n_ijk.Dot(n_ikl);
            sn = n_ikl.Cross(n_ijk).Dot(e_ik);
            if (sn != 0. || cs != 0.) {

              // Compute length of edge
              length = e_ik.Length();

              // Compute signed angle made up by normals in [-pi, pi]
              sn   /= length;
              angle = atan2(sn, cs);

              // Add mean curvature contribution of this edge (excl. factor 3)
              H = length * angle / A_sum;
              _Curvature->SetComponent(i, 0, _Curvature->GetComponent(i, 0) + H);
              _Curvature->SetComponent(k, 0, _Curvature->GetComponent(k, 0) + H);
              ++_Count[i], ++_Count[k];
            }
          }
        }
      }
    }
  }
};

// -----------------------------------------------------------------------------
/// Evaluate constraint penalty
struct Evaluate
{
  vtkDataArray *_MeanCurvature;
  double        _Penalty;

  Evaluate()
  :
    _Penalty(0.)
  {}

  Evaluate(const Evaluate &other, split)
  :
    _MeanCurvature(other._MeanCurvature),
    _Penalty(0.)
  {}

  void join(const Evaluate &other)
  {
    _Penalty += other._Penalty;
  }

  void operator ()(const blocked_range<int> &ptIds)
  {
    for (int ptId = ptIds.begin(); ptId != ptIds.end(); ++ptId) {
      _Penalty += abs(_MeanCurvature->GetComponent(ptId, 0));
    }
  }
};

// -----------------------------------------------------------------------------
/// Evaluate gradient of mean curvature
///
/// The following computes the actual gradient of the mean curvature as computed
/// by the vtkCurvatures::GetMeanCurvature function.
struct EvaluateGradient
{
  typedef MeanCurvatureConstraint::GradientType GradientType;

  vtkPolyData     *_Surface;
  vtkDataArray    *_Status;
  const EdgeTable *_EdgeTable;
  GradientType    *_Gradient;
  int             *_Count;

  /// Partial derivatives of cross product w.r.t. coordinates of a
  ///
  /// \note The Jacobian of the cross product w.r.t. the components of b
  ///       is equal to the negative Jacobian w.r.t. the components of a.
  ///
  /// \param[in] a First vector (unused).
  /// \param[in] b Second vector.
  ///
  /// \returns Partial derivatives of cross product w.r.t. components of vector a.
  static inline Matrix3x3 CrossJacobian(const Vector3 &a, const Vector3 &b)
  {
    return Matrix3x3(   0.,  b[2], -b[1],
                     -b[2],    0.,  b[0],
                      b[1], -b[0],    0.);
  }

  /// Compute gradient of vtkCurvatures::GetMeanCurvature
  void operator ()(const blocked_range<vtkIdType> &cellIds) const
  {
    vtkIdType i, j, k, l;
    vtkIdType numCellPts, *cellPts, cellPtIdx;
    vtkIdType numNborPts, *nborPts, nborPtIdx;
    double    v_i[3], v_j[3], v_k[3], v_l[3], cs2_plus_sn2, dangle_dsn, dangle_dcs;
    double    length, length2, A_ijk, A_ikl, A_sum, cs, sn, angle, H;
    Matrix3x3 dnJ_ijk, dnJ_ikl, dn_cross1, dn_cross2, dn_ijk, dn_ikl, T;
    Vector3   p_i, e_ij, e_ik, e_il, n_ijk, n_ikl, n_cross;
    Vector3   dA_sum, dlength, dcs, dsn, dangle, dH;

    vtkSmartPointer<vtkIdList> nborCellIds;
    nborCellIds = vtkSmartPointer<vtkIdList>::New();

    for (vtkIdType cellId = cellIds.begin(); cellId != cellIds.end(); ++cellId) {
      _Surface->GetCellPoints(cellId, numCellPts, cellPts);
      for (cellPtIdx = 0; cellPtIdx < numCellPts; ++cellPtIdx) {

        // Get start and end point indices of current cell edge
        i = cellPts[cellPtIdx];
        k = cellPts[(cellPtIdx + 1) % numCellPts];

        // Get neighboring face sharing this edge
        _Surface->GetCellEdgeNeighbors(cellId, i, k, nborCellIds);
        if (nborCellIds->GetNumberOfIds() == 1 && nborCellIds->GetId(0) > cellId) {

          // Get indices of two other points of this and the neighboring cell,
          // while preserving the order, i.e., orientation for both cells
          l = cellPts[(cellPtIdx + 2) % numCellPts];
          _Surface->GetCellPoints(nborCellIds->GetId(0), numNborPts, nborPts);
          if (numNborPts > 2) {
            nborPtIdx = 0;
            while (nborPts[nborPtIdx] != i) ++nborPtIdx;
            j = nborPts[(nborPtIdx + 1) % numNborPts];
            if (j == k) j = nborPts[(nborPtIdx + numNborPts - 1) % numNborPts];

            // Get vertex points
            _Surface->GetPoint(i, v_i);
            _Surface->GetPoint(j, v_j);
            _Surface->GetPoint(k, v_k);
            _Surface->GetPoint(l, v_l);

            // Compute required vector quantities from vertex points
            p_i   = Vector3(v_i);
            e_ij  = Vector3(v_j) - p_i;
            e_ik  = Vector3(v_k) - p_i;
            e_il  = Vector3(v_l) - p_i;
            n_ijk = e_ij.Cross(e_ik);
            n_ikl = e_ik.Cross(e_il);

            // Note: The area factor 1/2 is cancelled by the factor 2 of atan2
            A_ijk = n_ijk.Normalize();
            A_ikl = n_ikl.Normalize();
            A_sum = A_ijk + A_ikl;

            // Compute cross product of normal vectors (i.e., *after* Normalize)
            n_cross = n_ikl.Cross(n_ijk);

            // Compute cosine and sine of angle made up by the face normals
            cs = n_ijk.Dot(n_ikl);
            sn = n_cross.Dot(e_ik); // divided by l_ik inside if block
            if (sn != 0. || cs != 0.) {

              // Compute length of shared edge and its derivative
              length2 = e_ik.SquaredLength();
              length  = sqrt(length2);

              // Compute double angle using atan2 and the partial derivatives of atan2
              sn /= length;
              angle = atan2(sn, cs);

              cs2_plus_sn2 = cs * cs + sn * sn;
              dangle_dsn   =  cs / cs2_plus_sn2;
              dangle_dcs   = -sn / cs2_plus_sn2;

              //double z_norm     = sqrt(cs2_plus_sn2);
              //double inv_z_norm = 1. / z_norm;
              //double t          = (z_norm - cs) / sn;
              //double dt         = 2. / (1. + t * t);
              //dangle_dcs = - dt * (inv_z_norm + 1.) / sn;
              //dangle_dsn =   dt * (cs - inv_z_norm * sn - z_norm) / (sn * sn);

              // Compute other common terms
              dlength   = angle * e_ik / length; // w.r.t. v_k
              dn_cross1 =  CrossJacobian(n_ijk, n_ikl);
              dn_cross2 = -CrossJacobian(n_ikl, n_ijk);

              // Compute mean curvature (excl. factor 3)
              H = length * angle / A_sum;

              // Tensor product needed for partial derivatives of sine of angle
              T[0][0]           = e_ik[0] * e_ik[0] / length2 + 1.;
              T[0][1] = T[1][0] = e_ik[0] * e_ik[1] / length2;
              T[0][2] = T[2][0] = e_ik[0] * e_ik[2] / length2;
              T[1][1]           = e_ik[1] * e_ik[1] / length2 + 1.;
              T[1][2] = T[2][1] = e_ik[1] * e_ik[2] / length2;
              T[2][2]           = e_ik[2] * e_ik[2] / length2 + 1.;

              // Derivative of mean curvature w.r.t. v_i
              dnJ_ijk = Triangle::NormalDirectionJacobian(v_i, v_j, v_k);
              dnJ_ikl = Triangle::NormalDirectionJacobian(v_i, v_k, v_l);
              dn_ijk  = Triangle::NormalJacobian(n_ijk, dnJ_ijk);
              dn_ikl  = Triangle::NormalJacobian(n_ikl, dnJ_ikl);

              dA_sum  = n_ijk * dnJ_ijk;
              dA_sum += n_ikl * dnJ_ikl;
              dA_sum *= H;

              dsn  = e_ik * (dn_cross1 * dn_ijk + dn_cross2 * dn_ikl);
              dsn -= n_cross * T;
              dsn /= length;
              dsn *= dangle_dsn;

              dcs  = n_ikl * dn_ijk;
              dcs += n_ijk * dn_ikl;
              dcs *= dangle_dcs;

              dangle  = dsn;
              dangle += dcs;
              dangle *= length;

              dH  = -dlength;
              dH += dA_sum;
              dH += dangle;
              dH /= A_sum;

              _Gradient[i] += H * GradientType(dH);
              _Count   [i] += 1;

              // Derivative of mean curvature w.r.t. v_j
              dnJ_ijk = Triangle::NormalDirectionJacobian(v_j, v_k, v_i);
              dn_ijk  = Triangle::NormalJacobian(n_ijk, dnJ_ijk);

              dA_sum  = n_ijk * dnJ_ijk;
              dA_sum *= H;

              dsn  = e_ik * dn_cross1 * dn_ijk;
              dsn /= length;
              dsn *= dangle_dsn;

              dcs  = n_ikl * dn_ijk;
              dcs *= dangle_dcs;

              dangle  = dsn;
              dangle += dcs;
              dangle *= length;

              dH  = dA_sum;
              dH += dangle;
              dH /= A_sum;

              _Gradient[j] += H * GradientType(dH);
              _Count   [j] += 1;

              // Derivative of mean curvature w.r.t. v_k
              dnJ_ijk = Triangle::NormalDirectionJacobian(v_k, v_i, v_j);
              dnJ_ikl = Triangle::NormalDirectionJacobian(v_k, v_l, v_i);
              dn_ijk  = Triangle::NormalJacobian(n_ijk, dnJ_ijk);
              dn_ikl  = Triangle::NormalJacobian(n_ikl, dnJ_ikl);

              dA_sum  = n_ijk * dnJ_ijk;
              dA_sum += n_ikl * dnJ_ikl;
              dA_sum *= H;

              dsn  = e_ik * (dn_cross1 * dn_ijk + dn_cross2 * dn_ikl);
              dsn += n_cross * T;
              dsn /= length;
              dsn *= dangle_dsn;

              dcs  = n_ikl * dn_ijk;
              dcs += n_ijk * dn_ikl;
              dcs *= dangle_dcs;

              dangle  = dsn;
              dangle += dcs;
              dangle *= length;

              dH  = dlength;
              dH += dA_sum;
              dH += dangle;
              dH /= A_sum;

              _Gradient[k] += H * GradientType(dH);
              _Count   [k] += 1;

              // Derivative of mean curvature w.r.t. v_l
              dnJ_ikl = Triangle::NormalDirectionJacobian(v_l, v_i, v_k);
              dn_ikl  = Triangle::NormalJacobian(n_ikl, dnJ_ikl);

              dA_sum  = n_ikl * dnJ_ikl;
              dA_sum *= H;

              dsn  = e_ik * dn_cross2 * dn_ikl;
              dsn /= length;
              dsn *= dangle_dsn;

              dcs  = n_ijk * dn_ikl;
              dcs *= dangle_dcs;

              dangle  = dsn;
              dangle += dcs;
              dangle *= length;

              dH  = dA_sum;
              dH += dangle;
              dH /= A_sum;

              _Gradient[l] += H * GradientType(dH);
              _Count   [l] += 1;
            }
          }
        }
      }
    }
  }

//  /// Compute gradient of vtkCurvatures::GetMeanCurvature
//  void operator ()(const blocked_range<int> &ptIds) const
//  {
//    double     v_i[3], v_j[3], v_k[3], v_l[3];
//    double     l_ik, l2_ik, A_ijk, A_ikl, A_sum, cs_ik, sn_ik, a_ik, H_ik;
//    double     da_norm, da_ds, da_dc;
//    Matrix3x3  dnJ_ijk, dnJ_ikl, dn_cross1, dn_cross2, dn_ijk, dn_ikl, T;
//    Vector3    c, e_ij, e_ik, e_il, n_ijk, n_ikl, n_cross;
//    Vector3    dl_ik, dA_sum, dcs_ik, dsn_ik, da_ik, dH;
//    int        numAdjPts, adjPtIdx, i, j, k, l;
//    const int *adjPtIds;
//
//    for (i = ptIds.begin(); i != ptIds.end(); ++i) {
//      if (_Status && _Status->GetComponent(i, 0) == .0) continue;
//
//      _Surface->GetPoint(i, v_i);
//      _EdgeTable->GetAdjacentPoints(i, numAdjPts, adjPtIds);
//      if (numAdjPts > 0) {
//        for (adjPtIdx = 0; adjPtIdx < numAdjPts; ++adjPtIdx) {
//          k = adjPtIds[adjPtIdx];
//          if (i < k && GetEdgeNeighborPoints(_Surface, i, k, j, l) == 2 && j != -1 && l != -1) {
//
//            // Get points of adjacent triangle vertices
//            _Surface->GetPoint(j, v_j);
//            _Surface->GetPoint(k, v_k);
//            _Surface->GetPoint(l, v_l);
//
//            // Compute required vector quantities from vertex points
//            c[0] = v_i[0], c[1] = v_i[1], c[2] = v_i[2];
//            e_ij  = Vector3(v_j) - c;
//            e_ik  = Vector3(v_k) - c;
//            e_il  = Vector3(v_l) - c;
//            n_ijk = e_ij.Cross(e_ik);
//            n_ikl = e_ik.Cross(e_il);
//
//            // Note: The factor 1/2 is cancelled by the factor 2 of atan2
//            A_ijk = n_ijk.Normalize();
//            A_ikl = n_ikl.Normalize();
//            A_sum = A_ijk + A_ikl;
//            if (A_sum < 1e-6) continue;
//
//            // Compute cross product of normal vectors (i.e., *after* Normalize)
//            n_cross = n_ikl.Cross(n_ijk);
//
//            // Compute cosine and sine of angle made up by the face normals
//            cs_ik = n_ijk.Dot(n_ikl);
//            sn_ik = n_cross.Dot(e_ik); // divided by l_ik inside if block
//            if (abs(sn_ik) > 1e-3 && cs_ik > 1e-3) {
//
//              // Compute length of shared edge and its derivative
//              l2_ik = e_ik.SquaredLength();
//              l_ik  = sqrt(l2_ik);
//              if (l_ik < 1e-3) continue;
//
//              // Compute double angle using atan2 and the partial derivatives of atan2
//              sn_ik  /= l_ik;
//              a_ik    = atan2(sn_ik, cs_ik);
//              da_norm = cs_ik * cs_ik + sn_ik * sn_ik;
//              da_ds   =  cs_ik / da_norm;
//              da_dc   = -sn_ik / da_norm;
//
//              // Compute other common terms
//              dl_ik     = a_ik * e_ik / l_ik; // dl_ik / dv_k = - dl_ik / dv_i
//              dn_cross1 =  CrossJacobian(n_ijk, n_ikl);
//              dn_cross2 = -CrossJacobian(n_ikl, n_ijk);
//
//              // Compute mean curvature (excl. factor 3)
//              H_ik = l_ik * a_ik / A_sum;
//
//              // Tensor product needed for partial derivatives of sine of angle
//              T[0][0]           = e_ik[0] * e_ik[0] / l2_ik + 1.;
//              T[0][1] = T[1][0] = e_ik[0] * e_ik[1] / l2_ik;
//              T[0][2] = T[2][0] = e_ik[0] * e_ik[2] / l2_ik;
//              T[1][1]           = e_ik[1] * e_ik[1] / l2_ik + 1.;
//              T[1][2] = T[2][1] = e_ik[1] * e_ik[2] / l2_ik;
//              T[2][2]           = e_ik[2] * e_ik[2] / l2_ik + 1.;
//
//              // Derivative of mean curvature w.r.t. v_i
//              dnJ_ijk = Triangle::NormalDirectionJacobian(v_i, v_j, v_k);
//              dnJ_ikl = Triangle::NormalDirectionJacobian(v_i, v_k, v_l);
//              dn_ijk  = Triangle::NormalJacobian(n_ijk, dnJ_ijk);
//              dn_ikl  = Triangle::NormalJacobian(n_ikl, dnJ_ikl);
//
//              dA_sum  = n_ijk * dnJ_ijk;
//              dA_sum += n_ikl * dnJ_ikl;
//              dA_sum *= H_ik;
//
//              dsn_ik  = e_ik * (dn_cross1 * dn_ijk + dn_cross2 * dn_ikl);
//              dsn_ik -= n_cross * T;
//              dsn_ik /= l_ik;
//              dsn_ik *= da_ds;
//
//              dcs_ik  = n_ikl * dn_ijk;
//              dcs_ik += n_ijk * dn_ikl;
//              dcs_ik *= da_dc;
//
//              da_ik  = dsn_ik;
//              da_ik += dcs_ik;
//              da_ik *= l_ik;
//
//              dH  = -dl_ik;
//              dH += dA_sum;
//              dH += da_ik;
//              dH /= A_sum;
//
//              _Gradient[i] += H_ik * GradientType(dH);
//
//              // Derivative of mean curvature w.r.t. v_j
//              dnJ_ijk = Triangle::NormalDirectionJacobian(v_j, v_k, v_i);
//              dn_ijk  = Triangle::NormalJacobian(n_ijk, dnJ_ijk);
//
//              dA_sum  = n_ijk * dnJ_ijk;
//              dA_sum *= H_ik;
//
//              dsn_ik  = e_ik * dn_cross1 * dn_ijk;
//              dsn_ik /= l_ik;
//              dsn_ik *= da_ds;
//
//              dcs_ik  = n_ikl * dn_ijk;
//              dcs_ik *= da_dc;
//
//              da_ik  = dsn_ik;
//              da_ik += dcs_ik;
//              da_ik *= l_ik;
//
//              dH  = dA_sum;
//              dH += da_ik;
//              dH /= A_sum;
//
//              _Gradient[j] += H_ik * GradientType(dH);
//
//              // Derivative of mean curvature w.r.t. v_k
//              dnJ_ijk = Triangle::NormalDirectionJacobian(v_k, v_i, v_j);
//              dnJ_ikl = Triangle::NormalDirectionJacobian(v_k, v_l, v_i);
//              dn_ijk  = Triangle::NormalJacobian(n_ijk, dnJ_ijk);
//              dn_ikl  = Triangle::NormalJacobian(n_ikl, dnJ_ikl);
//
//              dA_sum  = n_ijk * dnJ_ijk;
//              dA_sum += n_ikl * dnJ_ikl;
//              dA_sum *= H_ik;
//
//              dsn_ik  = e_ik * (dn_cross1 * dn_ijk + dn_cross2 * dn_ikl);
//              dsn_ik += n_cross * T;
//              dsn_ik /= l_ik;
//              dsn_ik *= da_ds;
//
//              dcs_ik  = n_ikl * dn_ijk;
//              dcs_ik += n_ijk * dn_ikl;
//              dcs_ik *= da_dc;
//
//              da_ik  = dsn_ik;
//              da_ik += dcs_ik;
//              da_ik *= l_ik;
//
//              dH  = dl_ik;
//              dH += dA_sum;
//              dH += da_ik;
//              dH /= A_sum;
//
//              _Gradient[k] += H_ik * GradientType(dH);
//
//              // Derivative of mean curvature w.r.t. v_l
//              dnJ_ikl = Triangle::NormalDirectionJacobian(v_l, v_i, v_k);
//              dn_ikl  = Triangle::NormalJacobian(n_ikl, dnJ_ikl);
//
//              dA_sum  = n_ikl * dnJ_ikl;
//              dA_sum *= H_ik;
//
//              dsn_ik  = e_ik * dn_cross2 * dn_ikl;
//              dsn_ik /= l_ik;
//              dsn_ik *= da_ds;
//
//              dcs_ik  = n_ijk * dn_ikl;
//              dcs_ik *= da_dc;
//
//              da_ik  = dsn_ik;
//              da_ik += dcs_ik;
//              da_ik *= l_ik;
//
//              dH  = dA_sum;
//              dH += da_ik;
//              dH /= A_sum;
//
//              _Gradient[l] += H_ik * GradientType(dH);
//            }
//          }
//        }
//      }
//    }
//  }
};


#if USE_CURVATURE_WEIGHTED_SPRING_FORCE
// -----------------------------------------------------------------------------
/// Evaluate negative constraint force (i.e., gradient of constraint term)
///
/// Spring force designed to minimize curvature by weighting it based on mean curvature measures.
struct EvaluateWeightedSpringForce
{
  typedef MeanCurvatureConstraint::GradientType GradientType;

  vtkPoints       *_Points;
  vtkDataArray    *_Status;
  const EdgeTable *_EdgeTable;
  vtkDataArray    *_Normals;
  vtkDataArray    *_MeanCurvature;
  GradientType    *_Gradient;

  void operator ()(const blocked_range<int> &ptIds) const
  {
    double     c[3], p[3], d[3], f[3], m, H;
    int        numAdjPts;
    const int *adjPtIds;

    for (int ptId = ptIds.begin(); ptId != ptIds.end(); ++ptId) {
      if (_Status && _Status->GetComponent(ptId, 0) == .0) continue;
      _EdgeTable->GetAdjacentPoints(ptId, numAdjPts, adjPtIds);
      if (numAdjPts > 0) {
        // Magnitude of spring force based on mean curvature
        H = _MeanCurvature->GetComponent(ptId, 0);
//        if (H < 0.) m = 1. - SShapedMembershipFunction(-H, 0., .1);
//        else        m = SShapedMembershipFunction(H, 0., 1.);
        m = SShapedMembershipFunction(abs(H), 0., 1.);
        // Compute curvature weighted spring force
        _Points->GetPoint(ptId, c);
        f[0] = f[1] = f[2] = 0.;
        for (int i = 0; i < numAdjPts; ++i) {
          _Points->GetPoint(adjPtIds[i], p);
          vtkMath::Subtract(p, c, d);
          vtkMath::Add(f, d, f);
        }
        vtkMath::MultiplyScalar(f, 1. / numAdjPts);
        vtkMath::Normalize(f);
        _Gradient[ptId] = -m * GradientType(f);
      }
    }
  }
};
#endif // USE_CURVATURE_WEIGHTED_SPRING_FORCE


} // namespace MeanCurvatureConstraintUtils

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
void MeanCurvatureConstraint::CopyAttributes(const MeanCurvatureConstraint &other)
{
  _MaxMeanCurvature = other._MaxMeanCurvature;
}

// -----------------------------------------------------------------------------
MeanCurvatureConstraint::MeanCurvatureConstraint(const char *name, double weight)
:
  SurfaceConstraint(name, weight),
  _MaxMeanCurvature(NaN)
{
  _ParameterPrefix.push_back("Mean curvature ");
}

// -----------------------------------------------------------------------------
MeanCurvatureConstraint
::MeanCurvatureConstraint(const MeanCurvatureConstraint &other)
:
  SurfaceConstraint(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
MeanCurvatureConstraint &MeanCurvatureConstraint
::operator =(const MeanCurvatureConstraint &other)
{
  if (this != &other) {
    SurfaceConstraint::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
MeanCurvatureConstraint::~MeanCurvatureConstraint()
{
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
void MeanCurvatureConstraint::Initialize()
{
  // Initialize base class
  SurfaceConstraint::Initialize();

  // Add global (i.e., shared) point data array of computed surface curvatures
  #if USE_CURVATURE_WEIGHTED_SPRING_FORCE
    const bool global = true;
    AddPointData(SurfaceCurvature::MEAN, 1, VTK_FLOAT, global);
  #else // USE_CURVATURE_WEIGHTED_SPRING_FORCE
    AllocateCount(_NumberOfPoints);
  #endif // USE_CURVATURE_WEIGHTED_SPRING_FORCE
}

// -----------------------------------------------------------------------------
void MeanCurvatureConstraint::Update(bool gradient)
{
  // Update base class
  SurfaceConstraint::Update(gradient);

  if (_NumberOfPoints == 0) return;

  // Update mean curvature
  #if USE_CURVATURE_WEIGHTED_SPRING_FORCE
  vtkPolyData  * const surface        = DeformedSurface();
  vtkDataArray * const mean_curvature = PointData(SurfaceCurvature::MEAN);
  if (mean_curvature->GetMTime() < surface->GetMTime()) {

    SurfaceCurvature curv(SurfaceCurvature::Mean);
    curv.Input(surface);
    curv.EdgeTable(SharedEdgeTable());
    curv.VtkCurvaturesOn();
    curv.Run();

    MeshSmoothing smoother;
    smoother.Input(curv.Output());
    smoother.EdgeTable(SharedEdgeTable());
    smoother.SmoothPointsOff();
    smoother.SmoothArray(SurfaceCurvature::MEAN);
    smoother.NumberOfIterations(2);
    smoother.Run();

    vtkPointData * const smoothPD = smoother.Output()->GetPointData();
    mean_curvature->DeepCopy(smoothPD->GetArray(SurfaceCurvature::MEAN));
    mean_curvature->Modified();
  }
  #endif // USE_CURVATURE_WEIGHTED_SPRING_FORCE
}

// -----------------------------------------------------------------------------
double MeanCurvatureConstraint::Evaluate()
{
  #if USE_CURVATURE_WEIGHTED_SPRING_FORCE

    if (_NumberOfPoints == 0) return 0.;
    MeanCurvatureConstraintUtils::Evaluate eval;
    eval._MeanCurvature = PointData(SurfaceCurvature::MEAN);
    parallel_reduce(blocked_range<int>(0, _NumberOfPoints), eval);
    return eval._Penalty / _NumberOfPoints;

  #else // USE_CURVATURE_WEIGHTED_SPRING_FORCE

    const EdgeTable * const edgeTable = Edges();
    if (edgeTable->NumberOfEdges() == 0) return 0.;
    MeanCurvatureConstraintUtils::Evaluate eval;
    eval._Surface   = DeformedSurface();
    eval._Status    = Status();
    eval._EdgeTable = edgeTable;
    parallel_reduce(blocked_range<int>(0, _NumberOfPoints), eval);
    return 9. * eval._Penalty / edgeTable->NumberOfEdges();

  #endif // USE_CURVATURE_WEIGHTED_SPRING_FORCE
}

// -----------------------------------------------------------------------------
void MeanCurvatureConstraint
::EvaluateGradient(double *gradient, double step, double weight)
{
  if (_NumberOfPoints == 0 || _MaxMeanCurvature == 0.) return;

  memset(_Gradient, 0, _NumberOfPoints * sizeof(GradientType));

  #if USE_CURVATURE_WEIGHTED_SPRING_FORCE

    MeanCurvatureConstraintUtils::EvaluateWeightedSpringForce eval;
    eval._Points         = Points();
    eval._Status         = Status();
    eval._EdgeTable      = Edges();
    eval._Normals        = Normals();
    eval._MeanCurvature  = PointData(SurfaceCurvature::MEAN);
    eval._Gradient       = _Gradient;
    parallel_for(blocked_range<int>(0, _NumberOfPoints), eval);

  #else // USE_CURVATURE_WEIGHTED_SPRING_FORCE

    vtkPolyData * const surface = DeformedSurface();
    memset(_Count, 0, _NumberOfPoints * sizeof(int));
    MeanCurvatureConstraintUtils::EvaluateGradient eval;
    eval._Surface   = surface;
    eval._Status    = Status();
    eval._EdgeTable = Edges();
    eval._Gradient  = _Gradient;
    eval._Count     = _Count;
    eval(blocked_range<vtkIdType>(0, surface->GetNumberOfCells()));
    //parallel_for(blocked_range<int>(0, _NumberOfPoints), eval);
    for (int i = 0; i < _NumberOfPoints; ++i) {
      // - Factor 2 is from the derivative of the square function.
      // - Factor 3 is from the area divisor of the mean curvature.
      if (_Count[i] > 0) _Gradient[i] *= 6. / _Count[i];
    }

  #endif // USE_CURVATURE_WEIGHTED_SPRING_FORCE

  SurfaceConstraint::EvaluateGradient(gradient, step, weight / _NumberOfPoints);
}


} // namespace mirtk
