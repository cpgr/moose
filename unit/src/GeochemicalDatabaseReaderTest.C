//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef GEOCHEMICALDATABASEREADERTEST_H
#define GEOCHEMICALDATABASEREADERTEST_H

#include "gtest/gtest.h"

#include "GeochemicalDatabaseReader.h"

TEST(GeochemicalDatabaseReaderTest, reader)
{
  GeochemicalDatabaseReader database("data/geochem/sample.dat");

  // Set the primary, secondary and mineral species to be read
  std::vector<std::string> ps_names{"Ca++", "HCO3-", "H+"};
  std::vector<std::string> ss_names{"CO2(aq)", "CO3--", "CaCO3(aq)", "CaOH+", "OH-"};
  std::vector<std::string> ms_names{"Calcite"};

  database.setPrimarySpeciesNames(ps_names);
  database.setEquilibriumSpeciesNames(ss_names);
  database.setMineralSpeciesNames(ms_names);

  // Read the database and find these species
  database.read();

  // Get the temperature points from the database and compare with the expected
  // valules
  std::vector<Real> temperature_points;
  database.getTemperatures(temperature_points);

  std::vector<Real> temperature_points_gold{0.01, 25.0, 60.0, 100.0, 150.0, 200.0, 250.0, 300.0};
  EXPECT_EQ(temperature_points, temperature_points_gold);

  // Check the primary species
  std::vector<GeochemicalDatabasePrimarySpecies> ps;
  database.getPrimarySpecies(ps);

  bool unexpected_ps = false;
  for (auto s : ps)
  {
    if (s.name == "Ca++")
    {
      EXPECT_EQ(s.radius, 6);
      EXPECT_EQ(s.charge, 2);
      EXPECT_EQ(s.molar_mass, 40.078);
    }
    else if (s.name == "HCO3-")
    {
      EXPECT_EQ(s.radius, 4);
      EXPECT_EQ(s.charge, -1);
      EXPECT_EQ(s.molar_mass, 61.017);
    }
    else if (s.name == "H+")
    {
      EXPECT_EQ(s.radius, 9);
      EXPECT_EQ(s.charge, 1);
      EXPECT_EQ(s.molar_mass, 1.008);
    }
    else
      unexpected_ps = true;
  }

  // Check that no unexpected primary species were read
  EXPECT_FALSE(unexpected_ps);

  // Check the secondary species
  std::vector<GeochemicalDatabaseEquilibriumSpecies> ss;
  database.getEquilibriumSpecies(ss);

  bool unexpected_ss = false;
  for (auto s : ss)
  {
    if (s.name == "CO2(aq)")
    {
      std::vector<std::string> ps_gold{"H+", "H2O", "HCO3-"};
      std::vector<Real> sc_gold{1, -1, 1};
      std::vector<Real> logk_gold{
          -6.5804, -6.3447, -6.2684, -6.3882, -6.7235, -7.1969, -7.7868, -8.5280};
      EXPECT_EQ(s.debye_a, 3);
      EXPECT_EQ(s.charge, 0);
      EXPECT_EQ(s.molar_mass, 44.01);
      EXPECT_EQ(s.nspecies, 3);
      EXPECT_EQ(s.primary_species, ps_gold);
      EXPECT_EQ(s.stoichiometric_coeff, sc_gold);
      EXPECT_EQ(s.equilibrium_const, logk_gold);
    }
    else if (s.name == "CO3--")
    {
      std::vector<std::string> ps_gold{"H+", "HCO3-"};
      std::vector<Real> sc_gold{-1, 1};
      std::vector<Real> logk_gold{
          10.6241, 10.3288, 10.1304, 10.0836, 10.2003, 10.4648, 10.8707, 11.4638};
      EXPECT_EQ(s.debye_a, 4.5);
      EXPECT_EQ(s.charge, -2);
      EXPECT_EQ(s.molar_mass, 60.009);
      EXPECT_EQ(s.nspecies, 2);
      EXPECT_EQ(s.primary_species, ps_gold);
      EXPECT_EQ(s.stoichiometric_coeff, sc_gold);
      EXPECT_EQ(s.equilibrium_const, logk_gold);
    }
    else if (s.name == "CaCO3(aq)")
    {
      std::vector<std::string> ps_gold{"H+", "HCO3-", "Ca++"};
      std::vector<Real> sc_gold{-1, 1, 1};
      std::vector<Real> logk_gold{7.5021, 7.0017, 6.4516, 5.9636, 5.4683, 5.0185, 4.5355, 3.9118};
      EXPECT_EQ(s.debye_a, 3);
      EXPECT_EQ(s.charge, 0);
      EXPECT_EQ(s.molar_mass, 100.087);
      EXPECT_EQ(s.nspecies, 3);
      EXPECT_EQ(s.primary_species, ps_gold);
      EXPECT_EQ(s.stoichiometric_coeff, sc_gold);
      EXPECT_EQ(s.equilibrium_const, logk_gold);
    }
    else if (s.name == "CaOH+")
    {
      std::vector<std::string> ps_gold{"H2O", "Ca++", "H+"};
      std::vector<Real> sc_gold{1, 1, -1};
      std::vector<Real> logk_gold{500.0, 12.85, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0};
      EXPECT_EQ(s.debye_a, 4);
      EXPECT_EQ(s.charge, 1);
      EXPECT_EQ(s.molar_mass, 57.085);
      EXPECT_EQ(s.nspecies, 3);
      EXPECT_EQ(s.primary_species, ps_gold);
      EXPECT_EQ(s.stoichiometric_coeff, sc_gold);
      EXPECT_EQ(s.equilibrium_const, logk_gold);
    }
    else if (s.name == "OH-")
    {
      std::vector<std::string> ps_gold{"H2O", "H+"};
      std::vector<Real> sc_gold{1, -1};
      std::vector<Real> logk_gold{
          14.9398, 13.9951, 13.0272, 12.2551, 11.6308, 11.2836, 11.1675, 11.3002};
      EXPECT_EQ(s.debye_a, 3.5);
      EXPECT_EQ(s.charge, -1);
      EXPECT_EQ(s.molar_mass, 17.007);
      EXPECT_EQ(s.nspecies, 2);
      EXPECT_EQ(s.primary_species, ps_gold);
      EXPECT_EQ(s.stoichiometric_coeff, sc_gold);
      EXPECT_EQ(s.equilibrium_const, logk_gold);
    }
    else
      unexpected_ss = true;
  }

  // Check that no unexpected secondary species were read
  EXPECT_FALSE(unexpected_ss);

  // Check the mineral species
  std::vector<GeochemicalDatabaseMineralSpecies> ms;
  database.getMineralSpecies(ms);

  bool unexpected_ms = false;
  for (auto s : ms)
  {
    if (s.name == "Calcite")
    {
      std::vector<std::string> ps_gold{"H+", "HCO3-", "Ca++"};
      std::vector<Real> sc_gold{-1, 1, 1};
      std::vector<Real> logk_gold{
          2.2257, 1.8487, 1.3330, 0.7743, 0.0999, -0.5838, -1.3262, -2.2154};
      EXPECT_EQ(s.molar_volume, 36.934);
      EXPECT_EQ(s.molar_mass, 100.087);
      EXPECT_EQ(s.nspecies, 3);
      EXPECT_EQ(s.primary_species, ps_gold);
      EXPECT_EQ(s.stoichiometric_coeff, sc_gold);
      EXPECT_EQ(s.equilibrium_const, logk_gold);
    }
    else
      unexpected_ms = true;
  }

  // Check that no unexpected mineral species were read
  EXPECT_FALSE(unexpected_ms);
}

#endif // GEOCHEMICALDATABASEREADERTEST_H
