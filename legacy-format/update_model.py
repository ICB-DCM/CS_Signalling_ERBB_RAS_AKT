#!/usr/bin/env python3
import libsbml
import os
import amici


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


def check_sbml_consistency(sbml_document, check_units=False):
    """Check for SBML validity / consistency"""

    if not check_units:
        sbml_document.setConsistencyChecks(
            libsbml.LIBSBML_CAT_UNITS_CONSISTENCY, False)

    has_problems = sbml_document.checkConsistency()
    if has_problems:
        print_sbml_errors(sbml_document)
        print('WARNING: Generated invalid SBML model. Check messages above.')

    return has_problems


def print_sbml_errors(sbml_document,
                      minimum_severity=libsbml.LIBSBML_SEV_WARNING):
    """Print libsbml errors"""

    for error_idx in range(sbml_document.getNumErrors()):
        error = sbml_document.getError(error_idx)
        if error.getSeverity() >= minimum_severity:
            category = error.getCategoryAsString()
            severity = error.getSeverityAsString()
            message = error.getMessage()
            print(f'libSBML {severity} ({category}): {message}')


def globalize_parameters(sbml_model):
    """Turn all local parameters into global parameters"""
    for reaction in sbml_model.getListOfReactions():
        law = reaction.getKineticLaw()
        # copy first so we can delete in the following loop
        local_parameters = [local_parameter for local_parameter
                            in law.getListOfParameters()]
        for lp in local_parameters:
            # Create global parameter assuming names are globalized already
            p = sbml_model.createParameter()
            p.setId(lp.getId())
            p.setName(lp.getName())
            p.setConstant(lp.getConstant())
            p.setValue(lp.getValue())
            p.setUnits(lp.getUnits())

            # removeParameter, not removeLocalParameter!
            law.removeParameter(lp.getId())


def main():
    sbml_file_template = 'CS_Signalling_ERBB_RAS_AKT.xml'
    sbml_file_new = '../PEtab/CS_Signalling_ERBB_RAS_AKT_petab.xml'

    sbml_reader = libsbml.SBMLReader()
    sbml_document = sbml_reader.readSBML(sbml_file_template)
    sbml_model = sbml_document.getModel()

    print_sbml_errors(sbml_document)

    # fix units
    set_species_unit_to_nanomole(sbml_model)

    # Add observable and sigma
    add_observables(sbml_model)

    globalize_parameters(sbml_model)

    # Write updated model
    sbml_writer = libsbml.SBMLWriter()
    sbml_writer.writeSBMLToFile(sbml_document, sbml_file_new)

    # Load and check for errors
    sbml_reader = libsbml.SBMLReader()
    sbml_document = sbml_reader.readSBML(sbml_file_new)
    check_sbml_consistency(sbml_document)


if __name__ == '__main__':
    main()
