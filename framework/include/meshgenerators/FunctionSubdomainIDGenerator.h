//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MeshGenerator.h"

// Forward declarations
class FunctionSubdomainIDGenerator;
class GriddedData;

template <>
InputParameters validParams<FunctionSubdomainIDGenerator>();

/**
 * MeshGenerator for assigning a subdomain ID to all elements
 */
class FunctionSubdomainIDGenerator : public MeshGenerator
{
public:
  static InputParameters validParams();

  FunctionSubdomainIDGenerator(const InputParameters & parameters);

  std::unique_ptr<MeshBase> generate() override;

protected:
  /// the input mesh, which may be output by another mesh generator
  std::unique_ptr<MeshBase> & _input;

  /// direction where to look for value if interpolation order is constant
  MultiMooseEnum _direction;

  /// object to provide function evaluations at points on the grid
  std::unique_ptr<GriddedData> _gridded_data;
  /// dimension of the grid
  unsigned int _dim;

  /**
   * _axes specifies how to embed the grid into the MOOSE coordinate frame
   * if _axes[i] = 0 then the i_th axes of the grid lies along the MOOSE x direction
   * if _axes[i] = 1 then the i_th axes of the grid lies along the MOOSE y direction
   * if _axes[i] = 2 then the i_th axes of the grid lies along the MOOSE z direction
   * if _axes[i] = 3 then the i_th axes of the grid lies along the MOOSE time direction
   */
  std::vector<int> _axes;

  /// the grid
  std::vector<std::vector<Real>> _grid;

  /**
   * Operates on monotonically increasing in_arr.
   * Finds lower_x and upper_x which satisfy in_arr[lower_x] < x <= in_arr[upper_x].
   * End conditions: if x<in_arr[0] then lower_x = 0 = upper_x is returned
   *                 if x>in_arr[N-1] then lower_x = N-1 = upper_x is returned (N=size of in_arr)
   *
   * @param in_arr The monotonically increasing vector of real numbers
   * @param x The real value for which we want the neighbor indices
   * @param lower_x Upon return will contain lower_x specified above
   * @param upper_x Upon return will contain upper_x specified above
   */
  void getNeighborIndices(std::vector<Real> in_arr,
                          Real x,
                          unsigned int & lower_x,
                          unsigned int & upper_x) const;
};
