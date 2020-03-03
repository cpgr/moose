#!/usr/bin/env python
#* This file is part of the MOOSE framework
#* https://www.mooseframework.org
#*
#* All rights reserved, see COPYRIGHT for full restrictions
#* https://github.com/idaholab/moose/blob/master/COPYRIGHT
#*
#* Licensed under LGPL 2.1, please see LICENSE for details
#* https://www.gnu.org/licenses/lgpl-2.1.html

import os
import argparse
import numpy as np
import datetime

def command_line_options():
    """
    Command-line options for EQ3/6 thermodynamic reader and converter
    """
    parser = argparse.ArgumentParser(description="Utility to read data from EQ3/6 thermodynamic database and write to MOOSE thermodynamic database")
    parser.add_argument('-i', '--input', type=str, help="The EQ3/6 formatted input database", required=True)
    parser.add_argument('-o', '--output', type=str, default='eq36.dat', help="The name of the output file. Default filename is eq36.dat")
    return parser.parse_args()

def parseHeader(f):
    '''Reads the header of the EQ3/6 database'''
    dbheader = []
    f.seek(0)
    line = next(f)

    while not line.startswith('+--'):
        dbheader.append(line.strip())
        line = next(f)

    return dbheader

def parseMiscellaneous(f):
    '''Read the miscellaneous parameters of the EQ3/6 database'''

    line = next(f)
    temperatures = []
    dha = []
    dhb = []
    bdot = []

    while not line.startswith('+--'):
        # Read temperature points
        if line.startswith('temperatures'):
            line = next(f)
            while not line[0].isalpha():
                temperatures.append(line.split())
                line = next(f)

        # Read Debye-Huckel a
        if line.startswith('debye huckel a'):
            line = next(f)
            while not line[0].isalpha():
                dha.append(line.split())
                line = next(f)

        # Read Debye-Huckel b
        if line.startswith('debye huckel b'):
            line = next(f)
            while not line[0].isalpha():
                dhb.append(line.split())
                line = next(f)

        # Read Debye-Huckel bdot
        if line.startswith('bdot'):
            line = next(f)
            while not line[0].isalpha():
                bdot.append(line.split())
                line = next(f)

        line = next(f)

    tf = [item for sublist in temperatures for item in sublist]
    dhaf = [item for sublist in dha for item in sublist]
    dhbf = [item for sublist in dhb for item in sublist]
    bdotf = [item for sublist in bdot for item in sublist]

    return (tf, dhaf, dhbf, bdotf)

def parseBdot(f):
    '''Reads the bdot paramters of the EQ3/6 database'''
    a0 = {}
    line = next(f)

    while not line.startswith('+--'):
        a0[line.split()[0]] = line.split()[1]
        line = next(f)

    return a0

def parseBasisSpecies(f, species):
    '''Reads the basis species data of the EQ3/6 database'''

    # The current line is the last '+----' comment, so need to
    # go forward one line
    line = next(f)

    # Read all of the basis species until the auxiliary species
    while not line.startswith('auxiliary'):
        species_name = line.split()[0]
        species[species_name] = {}

        while not line.startswith('+--'):
            # Read molar weight
            if 'mol.wt.' in line:
                species[species_name]['molar_weight'] = line.split()[3]

            # Read a0
            if 'DHazero' in line:
                species[species_name]['a0'] = line.split()[3]

            # Read charge
            if 'charge' in line:
                species[species_name]['charge'] = line.split()[2]

            line = next(f)

        # Current line is '+---' after basis species, so proceed one line
        line = next(f)

    return

def parseAuxiliarySpecies(f, species, end_string):
    '''Reads the auxiliary species data of the EQ3/6 database'''
    line = next(f)

    # Read all of the basis species until end_string
    while not line.startswith(end_string):
        species_name = line.split()[0]
        species[species_name] = {}
        while not line.startswith('+--'):
            # Read molar weight
            if 'mol.wt.' in line:
                species[species_name]['molar_weight'] = line.split()[3]

            # Read a0
            if 'DHazero' in line:
                species[species_name]['a0'] = line.split()[3]

            # Read charge
            if 'charge' in line:
                species[species_name]['charge'] = line.split()[2]

            # Read constituent species
            if 'species in aqueous dissociation reaction' in line:
                reactants = []
                line = next(f)
                while not line.startswith('*'):
                    reactants.append(line.split())
                    line = next(f)

                reactants = [item for sublist in reactants for item in sublist]
                species[species_name]['elements'] = {}

                for i in range(0, len(reactants), 2):
                    species[species_name]['elements'][reactants[i+1]] = reactants[i]

            # Read equilibrium constant data
            if 'logK grid' in line:
                logk = []
                line = next(f)
                while not line.startswith('*'):
                    logk.append(line.split())
                    line = next(f)

                logk = [item for sublist in logk for item in sublist]
                species[species_name]['logk'] = logk

            line = next(f)

        # Current line is '+---' after basis species, so proceed one line
        line = next(f)

    return

def parseSolidSpecies(f, species):
    '''Reads the solid species data of the EQ3/6 database'''
    line = next(f)

    # Read all of the basis species until end_string
    while not line.startswith('liquids'):
        species_name = line.split()[0]
        species[species_name] = {}
        while not line.startswith('+--'):
            # Read molar weight
            if 'mol.wt.' in line:
                species[species_name]['molar_weight'] = line.split()[3]

            # Read molar volume
            if 'V0PrTr' in line:
                species[species_name]['molar_volume'] = line.split()[2]

            # Read constituent species
            if 'species in aqueous dissociation reaction' in line:
                reactants = []
                line = next(f)
                while not line.startswith('*'):
                    reactants.append(line.split())
                    line = next(f)

                reactants = [item for sublist in reactants for item in sublist]
                species[species_name]['elements'] = {}

                for i in range(0, len(reactants), 2):
                    species[species_name]['elements'][reactants[i+1]] = reactants[i]

            # Read equilibrium constant data
            if 'logK grid' in line:
                logk = []
                line = next(f)
                while not line.startswith('*'):
                    logk.append(line.split())
                    line = next(f)

                logk = [item for sublist in logk for item in sublist]
                species[species_name]['logk'] = logk

            line = next(f)

        line = next(f)

    return species


if __name__ == '__main__':

    # Command-line options
    opt = command_line_options()
    #opt.input = os.path.abspath(opt.input)
    opt.output

    # Open EQ3/6 database file for reading
    eq36f = open(opt.input, 'r')

    basis_species = {}
    secondary_species = {}
    solid_species = {}

    parsing_header = True

    for line in eq36f:

        # Parse header of EQ3/6 database
        if parsing_header:
            print "Parsing header"
            dbheader = parseHeader(eq36f)
            parsing_header = False

        # Parse the miscellaneous parameters data
        if line.startswith('Miscellaneous parameters') and next(eq36f).startswith('+--'):
            print "Parsing miscellaneous data"
            (temperatures, dha, dhb, bdot) = parseMiscellaneous(eq36f)

        # Parse the Debye-Huckel bdot data
        if line.startswith('bdot parameters') and next(eq36f).startswith('+--'):
            print "Parsing bdot data"
            a0 = parseBdot(eq36f)

        # Parse the primary basis species data
        if line.startswith('basis species') and next(eq36f).startswith('+--'):
            print "Parsing basis species"
            parseBasisSpecies(eq36f, basis_species)

            line = next(eq36f)

            # Parse the secondary basis species data
            #if line.startswith('auxiliary basis species') and next(eq36f).startswith('+--'):
            print "Parsing auxiliary basis species"
            parseAuxiliarySpecies(eq36f, secondary_species, 'aqueous')

            line = next(eq36f)

            # Parse the aqueous species data
            # if line.startswith('aqueous species') and next(eq36f).startswith('+--'):
            print "Parsing aqueous species"
            parseAuxiliarySpecies(eq36f, secondary_species, 'solids')
            # for key, value in aqueous_species.iteritems() :
            #     print key, value

            line = next(eq36f)

            # Parse the solid mineral species data
            #if line.startswith('solids') and next(eq36f).startswith('+--'):
            print "Parsing solid species"
            parseSolidSpecies(eq36f, solid_species)
            # for key, value in solid_species.iteritems() :
            #     print key, value

    # Close the EQ3/6 database filename
    eq36f.close()

    # Open output file for writing
    now = datetime.datetime.now()
    print "Finished parsing ", opt.input, "- writing data to ", opt.output

    outputf = open(opt.output, 'w')
    outputf.write('# Thermodynamic database for MOOSE chemical reactions module\n')
    outputf.write('# Created from EQ3/6 database ' + dbheader[0] + '\n')
    outputf.write('# Created ' + now.strftime("%H:%M %d-%m-%Y") + '\n')

    # Write the temperature points
    outputf.write('\n# Temperature points for log(K) in C\n')
    outputf.write('temperatures ')
    for i in range(0, len(temperatures)):
        outputf.write(temperatures[i] + ' ')

    # Write the Debye-Huckel parameters
    outputf.write('\n\n# Debye-Huckel activity coefficient parameters\n')
    outputf.write('Debye-Huckel parameters\n')
    outputf.write('a ')
    for i in range(0, len(dha)):
        outputf.write(dha[i] + ' ')

    outputf.write('\n')
    outputf.write('b ')
    for i in range(0, len(dhb)):
        outputf.write(dhb[i] + ' ')

    outputf.write('\n')
    outputf.write('bdot ')
    for i in range(0, len(bdot)):
        outputf.write(bdot[i] + ' ')

    # Write the primary species
    outputf.write('\n\n# Primary species: name, ion radius a0 (A), charge, molar weight (g/mol)\n')
    outputf.write('primary species\n')
    for item in sorted(basis_species):
        outputf.write(item + ' ')
        outputf.write(basis_species[item]['a0'] + ' ')
        outputf.write(basis_species[item]['charge'] + ' ')
        outputf.write(basis_species[item]['molar_weight'] + '\n')


    # Write the secondary species (auxiliary and aqueous)
    outputf.write('\n# Secondary species: name, ion radius a0 (A), charge, molar weight (g/mol), number of basis species, basis species and weights, log(K)\n')
    outputf.write('secondary species\n')
    for item in sorted(secondary_species):
        outputf.write(item + ' ')
        outputf.write(secondary_species[item]['a0'] + ' ')
        outputf.write(secondary_species[item]['charge'] + ' ')
        outputf.write(secondary_species[item]['molar_weight'] + ' ')

        # Write out the reactant species and stoichiometric coefficients
        num = len(secondary_species[item]['elements'].keys())
        outputf.write(str(num - 1) + ' ')
        for reactant in secondary_species[item]['elements']:
            if reactant != item:
                outputf.write(secondary_species[item]['elements'][reactant] + ' ')
                outputf.write(reactant + ' ')

        # Write the equilibrium constant values
        for i in range(0, len(secondary_species[item]['logk'])):
            outputf.write(secondary_species[item]['logk'][i] + ' ')

        outputf.write('\n')

    # Write the solid species
    outputf.write('\n# Solid mineral species: name, molar weight (g/mol), molar volume (cm^3/mol), number of basis species, basis species and weights, log(K)\n')
    outputf.write('minerals\n')
    for item in sorted(solid_species):
        outputf.write(item + ' ')
        outputf.write(solid_species[item]['molar_weight'] + ' ')
        outputf.write(solid_species[item]['molar_volume'] + ' ')

        # Write out the reactant species and stoichiometric coefficients
        num = len(solid_species[item]['elements'].keys())
        outputf.write(str(num - 1) + ' ')
        for reactant in solid_species[item]['elements']:
            if reactant != item:
                outputf.write(solid_species[item]['elements'][reactant] + ' ')
                outputf.write(reactant + ' ')

        # Write the equilibrium constant values
        for i in range(0, len(solid_species[item]['logk'])):
            outputf.write(solid_species[item]['logk'][i] + ' ')

        outputf.write('\n')

    # Close the output file
    outputf.close()
