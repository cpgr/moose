//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FunctionSubdomainIDGenerator.h"
#include "CastUniquePointer.h"
#include "GriddedData.h"

#include "libmesh/elem.h"

registerMooseObject("MooseApp", FunctionSubdomainIDGenerator);

defineLegacyParams(FunctionSubdomainIDGenerator);

InputParameters
FunctionSubdomainIDGenerator::validParams()
{
  InputParameters params = MeshGenerator::validParams();

  params.addRequiredParam<MeshGeneratorName>("input", "The mesh we want to modify");
  params.addParam<FileName>("data_file", "The data file of ids");
  MultiMooseEnum direction("left=0 right=1");
  params.addParam<MultiMooseEnum>(
      "direction", direction, "Direction to look to find value for each interpolation dimension.");

  params.addClassDescription("Sets all the elements of the input mesh to a unique subdomain ID.");

  return params;
}

FunctionSubdomainIDGenerator::FunctionSubdomainIDGenerator(const InputParameters & parameters)
  : MeshGenerator(parameters),
    _input(getMesh("input")),
    _direction(getParam<MultiMooseEnum>("direction")),
    _gridded_data(libmesh_make_unique<GriddedData>(getParam<FileName>("data_file"))),
    _dim(_gridded_data->getDim())
{
  _gridded_data->getAxes(_axes);
  _gridded_data->getGrid(_grid);

  // GriddedData does not require monotonicity of axes, but we do
  for (unsigned int i = 0; i < _dim; ++i)
    for (unsigned int j = 1; j < _grid[i].size(); ++j)
      if (_grid[i][j - 1] >= _grid[i][j])
        mooseError("FunctionSubdomainIDGenerator needs monotonically-increasing axis data.  Axis ",
                   i,
                   " contains non-monotonicity at value ",
                   _grid[i][j]);

  // GriddedData does not demand that each axis is independent, but we do
  std::set<int> s(_axes.begin(), _axes.end());
  if (s.size() != _dim)
    mooseError(
        "FunctionSubdomainIDGenerator needs the AXES to be independent.  Check the AXIS lines in "
        "your data file.");
}

std::unique_ptr<MeshBase>
FunctionSubdomainIDGenerator::generate()
{
  std::unique_ptr<MeshBase> mesh = std::move(_input);

  for (auto & elem : mesh->element_ptr_range())
  {
    std::vector<unsigned int> left(_dim);
    std::vector<unsigned int> right(_dim);
    std::vector<unsigned int> arg(_dim);
    Point pt = elem->vertex_average();

    for (unsigned int i = 0; i < _dim; ++i)
    {
      getNeighborIndices(_grid[i], pt(i), left[i], right[i]);
      if (_direction.get(i) == 0)
        arg[i] = left[i];
      else
        arg[i] = right[i];
    }

    elem->subdomain_id() = _gridded_data->evaluateFcn(arg);
  }

  return dynamic_pointer_cast<MeshBase>(mesh);
}

void
FunctionSubdomainIDGenerator::getNeighborIndices(std::vector<Real> in_arr,
                                                 Real x,
                                                 unsigned int & lower_x,
                                                 unsigned int & upper_x) const
{
  int N = in_arr.size();
  if (x <= in_arr[0])
  {
    lower_x = 0;
    upper_x = 0;
  }
  else if (x >= in_arr[N - 1])
  {
    lower_x = N - 1;
    upper_x = N - 1;
  }
  else
  {
    // returns up which points at the first element in inArr that is not less than x
    std::vector<double>::iterator up = std::lower_bound(in_arr.begin(), in_arr.end(), x);

    // std::distance returns std::difference_type, which can be negative in theory, but
    // in this context will always be >=0.  Therefore the explicit cast is just to shut
    // the compiler up.
    upper_x = static_cast<unsigned int>(std::distance(in_arr.begin(), up));
    if (in_arr[upper_x] == x)
      lower_x = upper_x;
    else
      lower_x = upper_x - 1;
  }
}
