//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "PorousFlowDiracKernelScalarBase.h"
#include "ReporterInterface.h"

/**
 * Approximates a borehole by a sequence of Dirac Points
 *
 * The base is PorousFlowDiracKernelScalarBase rather than DiracKernel so that a line sink can hold
 * a bottom-hole pressure that is itself an unknown.  Without a "scalar_variable" that base is a
 * plain DiracKernel, so line sinks that do not use one are unaffected.
 */
class PorousFlowLineGeometry : public PorousFlowDiracKernelScalarBase, public ReporterInterface
{
public:
  /**
   * Creates a new PorousFlowLineGeometry
   * This reads the file containing the lines of the form
   * weight x y z
   * that defines the line geometry.
   * It also calculates segment-lengths between the points
   */
  static InputParameters validParams();

  PorousFlowLineGeometry(const InputParameters & parameters);

  void virtual initialSetup() override;

protected:
  /// Line length.  This is only used if there is only one borehole point
  const Real _line_length;

  /// Line direction.  This is only used if there is only one borehole point
  const RealVectorValue _line_direction;

  /**
   * File defining the geometry of the borehole.   Each row has format
   * weight x y z
   * and the list of such points defines a polyline that is the line sink
   */
  const std::string _point_file;

  /// Radii of the borehole
  std::vector<Real> _rs;

  /// x points of the borehole
  std::vector<Real> _xs;

  /// y points of the borehole
  std::vector<Real> _ys;

  /// z points of borehole
  std::vector<Real> _zs;

  /// The bottom point of the borehole (where bottom_pressure is defined)
  Point _bottom_point;

  /// 0.5*(length of polyline segments between points)
  std::vector<Real> _half_seg_len;

  /// Length of the entire polyline, which is completionLength() summed over every point
  Real _total_line_length;

  /**
   * The length of line attributed to one Dirac point: half the distance back to the previous
   * point plus half the distance forward to the next, so that summing over all the points
   * recovers the length of the line
   * @param dirac_ptid The ID assigned to the point by addPoints()
   */
  Real completionLength(unsigned dirac_ptid) const;

  /**
   * The fraction of the line attributed to one Dirac point, so that the fractions of all the
   * points sum to unity.  This is the share of a prescribed total, such as an injection rate,
   * that this point should carry.  A line with a single completion carries all of it whatever
   * line_length says: unlike the Peaceman well constant, a prescribed rate does not scale with
   * the length of the completion, so line_length may be left at its default here.
   * @param dirac_ptid The ID assigned to the point by addPoints()
   */
  Real completionFraction(unsigned dirac_ptid) const;

  /// Add Dirac Points to the line sink
  virtual void addPoints() override;

  /// regenerate points in each cell if using line_base
  virtual void meshChanged() override;

  /// Reads a space-separated line of floats from ifs and puts in myvec
  bool parseNextLineReals(std::ifstream & ifs, std::vector<Real> & myvec);

  ///@{
  /// reporter input alternative (to the point file and line_base data)
  const std::vector<Real> * const _x_coord;
  const std::vector<Real> * const _y_coord;
  const std::vector<Real> * const _z_coord;
  const std::vector<Real> * const _weight;
  const bool _usingReporter;
  ///@}
private:
  void calcLineLengths();
  void regenPoints();

  /// alternative (to the point file data) line weight and start point.
  std::vector<Real> _line_base;
};
