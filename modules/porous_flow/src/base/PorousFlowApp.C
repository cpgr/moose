//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PorousFlowApp.h"
#include "Moose.h"
#include "SolidMechanicsApp.h"
#include "AppFactory.h"
#include "MooseSyntax.h"
#include "FluidPropertiesApp.h"
#include "ChemicalReactionsApp.h"

InputParameters
PorousFlowApp::validParams()
{
  InputParameters params = MooseApp::validParams();

  params.set<bool>("automatic_automatic_scaling") = false;
  params.set<bool>("use_legacy_material_output") = false;
  params.set<bool>("use_legacy_initial_residual_evaluation_behavior") = false;

  return params;
}

registerKnownLabel("PorousFlowApp");

PorousFlowApp::PorousFlowApp(const InputParameters & parameters) : MooseApp(parameters)
{
  PorousFlowApp::registerAll(_factory, _action_factory, _syntax);
}

PorousFlowApp::~PorousFlowApp() {}

void
PorousFlowApp::registerApps()
{
  registerApp(PorousFlowApp);

  SolidMechanicsApp::registerApps();
  FluidPropertiesApp::registerApps();
  ChemicalReactionsApp::registerApps();
}

void
PorousFlowApp::registerAll(Factory & f, ActionFactory & af, Syntax & syntax)
{
  SolidMechanicsApp::registerAll(f, af, syntax);
  FluidPropertiesApp::registerAll(f, af, syntax);
  ChemicalReactionsApp::registerAll(f, af, syntax);
  Registry::registerObjectsTo(f, {"PorousFlowApp"});
  Registry::registerActionsTo(af, {"PorousFlowApp"});

  registerSyntax("CoupledInjectionProductionPhysics", "Physics/CoupledInjectionProduction/*");

  registerSyntaxTask("PorousFlowUnsaturated", "PorousFlowUnsaturated", "add_user_object");
  registerSyntaxTask("PorousFlowUnsaturated", "PorousFlowUnsaturated", "add_kernel");
  registerSyntaxTask("PorousFlowUnsaturated", "PorousFlowUnsaturated", "add_material");
  registerSyntaxTask("PorousFlowUnsaturated", "PorousFlowUnsaturated", "add_aux_variable");
  registerSyntaxTask("PorousFlowUnsaturated", "PorousFlowUnsaturated", "add_aux_kernel");

  registerSyntaxTask("PorousFlowFullySaturated", "PorousFlowFullySaturated", "add_user_object");
  registerSyntaxTask("PorousFlowFullySaturated", "PorousFlowFullySaturated", "add_kernel");
  registerSyntaxTask("PorousFlowFullySaturated", "PorousFlowFullySaturated", "add_material");
  registerSyntaxTask("PorousFlowFullySaturated", "PorousFlowFullySaturated", "add_aux_variable");
  registerSyntaxTask("PorousFlowFullySaturated", "PorousFlowFullySaturated", "add_aux_kernel");

  registerSyntaxTask("PorousFlowBasicTHM", "PorousFlowBasicTHM", "add_user_object");
  registerSyntaxTask("PorousFlowBasicTHM", "PorousFlowBasicTHM", "add_kernel");
  registerSyntaxTask("PorousFlowBasicTHM", "PorousFlowBasicTHM", "add_material");
  registerSyntaxTask("PorousFlowBasicTHM", "PorousFlowBasicTHM", "add_aux_variable");
  registerSyntaxTask("PorousFlowBasicTHM", "PorousFlowBasicTHM", "add_aux_kernel");

  syntax.registerActionSyntax("PorousFlowAddBCAction", "Modules/PorousFlow/BCs/*");

  registerMooseObjectTask("add_porous_flow_bc", PorousFlowSinkBC, false);
  addTaskDependency("add_porous_flow_bc", "add_bc");
  addTaskDependency("resolve_optional_materials", "add_porous_flow_bc");

  // Wells are operational forcing on the PorousFlow equations rather than PDE terms, so [Wells] is
  // a top-level block that reads only the PorousFlowDictator: it works identically with
  // Physics-built and hand-built inputs
  syntax.registerActionSyntax("PorousFlowWellAction", "Wells/*", "add_porous_flow_well_variable");
  syntax.registerActionSyntax("PorousFlowWellAction", "Wells/*", "add_porous_flow_well_gate");
  syntax.registerActionSyntax("PorousFlowWellAction", "Wells/*", "add_porous_flow_well");

  // A well whose bottom-hole pressure is solved for adds a SCALAR variable to the nonlinear system,
  // which has to happen before the variable window closes.  The gate task below is too late for
  // that: it sits after add_aux_variable, which itself follows add_variable
  registerMooseObjectTask("add_porous_flow_well_variable", PorousFlowWellBase, false);
  addTaskDependency("add_variable", "add_porous_flow_well_variable");
  addTaskDependency("add_scalar_kernel", "add_porous_flow_well_variable");

  // A gated production well adds an AuxVariable, which must appear in the normal variable window:
  // early enough for the relationship managers to see it and for its degrees of freedom to be
  // distributed with everything else
  registerMooseObjectTask("add_porous_flow_well_gate", PorousFlowWellBase, false);
  addTaskDependency("add_porous_flow_well_gate", "add_aux_variable");
  addTaskDependency("add_aux_kernel", "add_porous_flow_well_gate");

  registerMooseObjectTask("add_porous_flow_well", PorousFlowWellBase, false);
  addTaskDependency("add_porous_flow_well", "add_porous_flow_well_gate");
  addTaskDependency("add_porous_flow_well", "add_dirac_kernel");
  addTaskDependency("resolve_optional_materials", "add_porous_flow_well");

  // Task dependency and syntax for action to automatically add PorousFlow materials
  registerSyntax("PorousFlowAddMaterialAction", "Materials");

  // Task dependency and syntax for action to automatically add PorousFlowJoiner materials
  registerTask("add_joiners", /*is_required=*/false);
  addTaskDependency("add_joiners", "add_material");
  addTaskDependency("add_master_action_material", "add_joiners");

  registerSyntaxTask("PorousFlowAddMaterialJoiner", "Materials", "add_joiners");
}

extern "C" void
PorousFlowApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  PorousFlowApp::registerAll(f, af, s);
}
extern "C" void
PorousFlowApp__registerApps()
{
  PorousFlowApp::registerApps();
}
