#!/usr/bin/env python3
import libsbml
import os
import amici
import petab
import re

def set_species_unit_to_nanomole(sbml_model):
    # create nanomole
    unitDefinition = sbml_model.createUnitDefinition()
    unitDefinition.setId('nmol')
    unitDefinition.setId('nano mole')
    unit = unitDefinition.createUnit()
    unit.setKind(libsbml.UNIT_KIND_MOLE)
    unit.setExponent(1)
    unit.setScale(-9)
    unit.setMultiplier(1)

    for s in sbml_model.getListOfSpecies():
        s.setUnits('nmol')


def create_parameter(model, id, constant, value, units):
    k = model.createParameter()
    k.setId(id)
    k.setName(id)
    k.setConstant(constant)
    k.setValue(value)
    k.setUnits(units)
    return k


def create_assigment_rule(model, name, formula):
    rule = model.createAssignmentRule()
    rule.setId(name)
    rule.setName(name)
    rule.setVariable(name)
    rule.setFormula(formula)
    return rule


def add_sigma(model, observable):
    p = create_parameter(model, f'sigma_{observable}', False, 1.0,
                         'dimensionless')
    p = create_parameter(model, f'noiseParameter1_{observable}', True, 1.0,
                         'dimensionless')
    rule = create_assigment_rule(model, f'sigma_{observable}',
                                 f'noiseParameter1_{observable}')


def add_observables(sbml_model):
    # Remove unused protein observables
    par_list = [p.getId() for p in sbml_model.getListOfParameters()]
    for p in par_list:
        if p.startswith('offset_protein_common_'):
            obs = 'Mclp_' + p[len('offset_protein_common_'):]
            pp = sbml_model.getParameter(p)
            pp.setId(f'observableParameter1_{obs}')
            pp.setValue(0.0)
            create_parameter(sbml_model, f'observableParameter2_{obs}', False,
                             1.0, 'dimensionless')
            r = sbml_model.getRuleByVariable(f'observable_{obs}')
            r.setFormula(re.sub(p,
                                f'observableParameter1_{obs}',
                                r.getFormula()))
            r.setFormula(re.sub('offset_protein_total_genotypespecific',
                                f'observableParameter2_{obs}',
                                r.getFormula()))

            add_sigma(sbml_model, obs)

        if p.endswith('_sigma'):
            sbml_model.removeParameter(p)

        if p.startswith('weight_proliferation_') \
                or p == 'scaling_proliferation_reference':
            sbml_model.getParameter(p).setValue(1.0)

        if p.startswith('weight_proliferation_') \
                or p == 'scaling_proliferation_reference':
            sbml_model.getParameter(p).setValue(1.0)

    sbml_model.removeParameter('protein_sigma')
    sbml_model.removeParameter('proliferation_sigma')
    sbml_model.removeParameter('offset_protein_total_genotypespecific')

    sbml_model.getRuleByVariable('observable_proliferation').setVariable('observable_viability')
    sbml_model.getParameter('observable_proliferation').setId('observable_viability')
    sbml_model.getParameter('scaling_proliferation_reference').setId(
        'observableParameter1_viability')
    add_sigma(sbml_model, 'viability')
    r = sbml_model.getRuleByVariable('observable_viability')
    r.setFormula(re.sub('scaling_proliferation_reference',
                        'observableParameter1_viability', r.getFormula()))


def gene_specific_scaling_non_const(sbml_model):
    """Make *_k_GeneSpecificScaling parameter non-const"""
    for p in sbml_model.getListOfParameters():
        if p.getId().endswith('_k_GeneSpecificScaling'):
            p.setConstant(False)


def main():
    sbml_file_template = 'ERBB_RAS_AKT_Drugs_r389936_withObs_withSigmaPerProtein.sbml'
    sbml_file_new = 'CS_Signalling_ERBB_RAS_AKT_CCLE_MCLP_petab.xml'

    sbml_reader = libsbml.SBMLReader()
    sbml_document = sbml_reader.readSBML(sbml_file_template)
    sbml_model = sbml_document.getModel()

    petab.log_sbml_errors(sbml_document)

    # fix units
    set_species_unit_to_nanomole(sbml_model)

    # Add observable and sigma
    add_observables(sbml_model)

    petab.globalize_parameters(sbml_model)

    gene_specific_scaling_non_const(sbml_model)

    # Write updated model
    sbml_writer = libsbml.SBMLWriter()
    sbml_writer.writeSBMLToFile(sbml_document, sbml_file_new)

    # Load and check for errors
    sbml_reader = libsbml.SBMLReader()
    sbml_document = sbml_reader.readSBML(sbml_file_new)
    petab.is_sbml_consistent(sbml_document)


if __name__ == '__main__':
    main()
