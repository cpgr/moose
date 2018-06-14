//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef POROUSFLOWBRINEMETHANETEST_H
#define POROUSFLOWBRINEMETHANETEST_H

#include "gtest_include.h"

#include "FEProblem.h"
#include "AppFactory.h"
#include "GeneratedMesh.h"
#include "PorousFlowCapillaryPressureVG.h"
#include "PorousFlowBrineMethane.h"
#include "BrineFluidProperties.h"
#include "Water97FluidProperties.h"
#include "NaClFluidProperties.h"
#include "MethaneFluidProperties.h"
#include "MooseApp.h"
#include "Utils.h"

class PorousFlowBrineMethaneTest : public ::testing::Test
{
protected:
  void SetUp()
  {
    const char * argv[] = {"foo", NULL};

    _app = AppFactory::createAppShared("MooseUnitApp", 1, (char **)argv);
    _factory = &_app->getFactory();

    registerObjects(*_factory);
    buildObjects();
  }

  void registerObjects(Factory & factory)
  {
    registerUserObject(PorousFlowCapillaryPressureVG);
    registerUserObject(BrineFluidProperties);
    registerUserObject(Water97FluidProperties);
    registerUserObject(NaClFluidProperties);
    registerUserObject(MethaneFluidProperties);
    registerUserObject(PorousFlowBrineMethane);
  }

  void buildObjects()
  {
    InputParameters mesh_params = _factory->getValidParams("GeneratedMesh");
    mesh_params.set<MooseEnum>("dim") = "3";
    mesh_params.set<std::string>("name") = "mesh";
    mesh_params.set<std::string>("_object_name") = "name1";
    _mesh = libmesh_make_unique<GeneratedMesh>(mesh_params);

    InputParameters problem_params = _factory->getValidParams("FEProblem");
    problem_params.set<MooseMesh *>("mesh") = _mesh.get();
    problem_params.set<std::string>("name") = "problem";
    problem_params.set<std::string>("_object_name") = "name2";
    _fe_problem = libmesh_make_unique<FEProblem>(problem_params);

    InputParameters pc_params = _factory->getValidParams("PorousFlowCapillaryPressureVG");
    pc_params.set<Real>("m") = 0.5;
    pc_params.set<Real>("alpha") = 0.1;
    _fe_problem->addUserObject("PorousFlowCapillaryPressureVG", "pc", pc_params);
    _pc = &_fe_problem->getUserObject<PorousFlowCapillaryPressureVG>("pc");

    InputParameters brine_params = _factory->getValidParams("BrineFluidProperties");
    _fe_problem->addUserObject("BrineFluidProperties", "brine_fp", brine_params);
    _brine_fp = &_fe_problem->getUserObject<BrineFluidProperties>("brine_fp");

    InputParameters water_params = _factory->getValidParams("Water97FluidProperties");
    _fe_problem->addUserObject("Water97FluidProperties", "water_fp", water_params);
    _water_fp = &_fe_problem->getUserObject<Water97FluidProperties>("water_fp");

    InputParameters methane_params = _factory->getValidParams("MethaneFluidProperties");
    _fe_problem->addUserObject("MethaneFluidProperties", "methane_fp", methane_params);
    _methane_fp = &_fe_problem->getUserObject<MethaneFluidProperties>("methane_fp");

    InputParameters uo_params = _factory->getValidParams("PorousFlowBrineMethane");
    uo_params.set<UserObjectName>("brine_fp") = "brine_fp";
    uo_params.set<UserObjectName>("methane_fp") = "methane_fp";
    uo_params.set<UserObjectName>("capillary_pressure") = "pc";
    _fe_problem->addUserObject("PorousFlowBrineMethane", "fp", uo_params);
    _fp = &_fe_problem->getUserObject<PorousFlowBrineMethane>("fp");
  }

  Real celsiusToKelvin(Real c) {return c + 273.15;}

  Real barToPascal(Real b) {return b * 1.0e5;}

  MooseAppPtr _app;
  std::unique_ptr<MooseMesh> _mesh;
  std::unique_ptr<FEProblem> _fe_problem;
  Factory * _factory;
  const PorousFlowCapillaryPressureVG * _pc;
  const PorousFlowBrineMethane * _fp;
  const BrineFluidProperties * _brine_fp;
  const Water97FluidProperties * _water_fp;
  const MethaneFluidProperties * _methane_fp;

  // Perturbation magnitude for the computation of derivatives
  const Real dp = 1.0e-2;  // Pressure in Pascal
  const Real dT = 1.0e-6;  // Temperature in Kelvin
};

#endif // POROUSFLOWBRINEMETHANETEST_H
