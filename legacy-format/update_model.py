#!/usr/bin/env python3
import libsbml
import os
import amici
import petab


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
    # Remove old proliferation observable
    sbml_model.removeRuleByVariable('proliferation_score')
    sbml_model.removeParameter('proliferation_score')
    for i in range(16):
        sbml_model.removeParameter(f'k{i}')


    # Add proliferation observable and sigma
    observable_formula = """observableParameter1_proliferation * 
    ( SP_64_6 * weight_proliferation_SP_64_6 
    + SP_258_6 * weight_proliferation_SP_258_6 
    + SP_259_6 * weight_proliferation_SP_259_6 
    + SP_382_6 * weight_proliferation_SP_382_6 
    + SP_477_6 * weight_proliferation_SP_477_6 
    + SP_478_6 * weight_proliferation_SP_478_6 
    + SP_1512_6 * weight_proliferation_SP_1512_6 
    + SP_1676_6 * weight_proliferation_SP_1676_6 
    + SP_18461_6 * weight_proliferation_SP_18461_6 
    + SP_22030_6 * weight_proliferation_SP_22030_6 
    + SP_22032_6 * weight_proliferation_SP_22032_6 
    + SP_22729_6 * weight_proliferation_SP_22729_6 ) 
    / ( SP_493_6 * weight_proliferation_SP_493_6 
    + SP_494_6 * weight_proliferation_SP_494_6 
    + SP_495_6 * weight_proliferation_SP_495_6 
    + SP_496_6 * weight_proliferation_SP_496_6 + 1 )"""
    observable_formula = observable_formula.replace('\n', ' ')

    for s in observable_formula.split(' '):
        if s.startswith('weight_'):
            create_parameter(sbml_model, s, False, 1.0, 'dimensionless')
    create_parameter(sbml_model, 'observableParameter1_proliferation', False,
                     1.0, 'dimensionless')

    observable_id = 'observable_proliferation'
    p = create_parameter(sbml_model, observable_id, False, 1.0, 'dimensionless')
    rule = create_assigment_rule(sbml_model, observable_id, observable_formula)
    add_sigma(sbml_model, f'proliferation')


def gene_specific_scaling_non_const(sbml_model):
    """Make *_k_GeneSpecificScaling parameter non-const"""
    for p in sbml_model.getListOfParameters():
        if p.getId().endswith('_k_GeneSpecificScaling'):
            p.setConstant(False)

def main():
    sbml_file_template = 'CS_Signalling_ERBB_RAS_AKT.xml'
    sbml_file_new = '../PEtab/CS_Signalling_ERBB_RAS_AKT_petab.xml'

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
