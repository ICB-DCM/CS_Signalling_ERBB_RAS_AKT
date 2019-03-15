#!/usr/bin/env python3
import petab
import pandas as pd
import re
import numpy as np
import libsbml


def update_parameter_name(old_name):
    # p_r_23811_k_RPKM2protein_reaction_23811
    # -> reaction_23811_r23811_k_RPKM2protein
    # if no match, return as is
    # new_name = re.sub(r'p_r_(\d+)_(.*)_(reaction_\d+)',
    # r'\3_r\1_\2', old_name)
    #new_name = re.sub(r'p_r_(\d+_.*)_(reaction_\d+)', r'r\1', old_name)
    new_name = re.sub(r'(p_r_\d+_.*)_(reaction_\d+)', r'\1', old_name)
    return new_name


def update_condition_table(condition_file_template, condition_file, sbml_model):
    template_df = pd.read_csv(condition_file_template, sep='\t')
    template_df.ID = template_df.ID.apply(update_parameter_name)
    template_df.set_index(['ID'], inplace=True)
    template_df = template_df.transpose()
    template_df.index.name = 'conditionId'

    # Remove extra states. Only keep the ones marked constant in the model
    for parameter_id in template_df:
        if parameter_id.endswith('_k_GeneSpecificScaling'):
            del template_df[parameter_id]
            continue

        if not parameter_id.startswith('SP_'):
            continue

        state_id = parameter_id
        species = sbml_model.getSpecies(state_id)
        if species and not species.getConstant():
            del template_df[state_id]

    template_df.to_csv(condition_file, sep='\t', index=True)


def update_measurement_table(measurement_file_template, measurement_file):
    measurement_df = pd.read_csv(measurement_file_template, sep='\t')
    del measurement_df['conditionRef']
    # template_df.set_index(['conditionId'], inplace=True)
    measurement_df.rename(columns={'observable': 'observableId',
                                   'condition': 'simulationConditionId',
                                   'sigma': 'noiseParameters',
                                   'scalingParameter': 'observableParameters'},
                          inplace=True)
    # shorten
    # scaling_proliferation_reference_genotypespecific_TUMOR-A2058-cellline-01-01
    measurement_df.observableParameters = \
        measurement_df.observableParameters.apply(
            lambda x: x.replace('reference_genotypespecific_', '')
        )

    measurement_df.observableId = \
        measurement_df.observableId.apply(
            lambda x: re.sub('proliferation',
                             'viability', x)
        )

    measurement_df.observableId = \
        measurement_df.observableId.apply(
            lambda x: re.sub('proliferation',
                             'viability', x)
        )

    for i, row in measurement_df.iterrows():
        spl = row.observableParameters.split(',')
        if len(spl) > 1:
            # two-dataset table
            obs_pars = []
            for p in spl:
                if p.find('sigma') >= 0:
                    measurement_df.loc[i, 'noiseParameters'] = p
                else:
                    obs_pars.append(p)
            measurement_df.loc[i, 'observableParameters'] = ';'.join(obs_pars)

    measurement_df['preequilibrationConditionId'] = np.nan
    measurement_df['observableTransformation'] = np.nan
    measurement_df.to_csv(measurement_file, sep='\t', index=False)


def create_parameter_table(sbml_file, condition_file, measurement_file,
                           parameter_file):
    """Create PEtab parameter table"""

    problem = petab.Problem.from_files(
        sbml_file=sbml_file,
        condition_file=condition_file,
        measurement_file=measurement_file)

    df = problem.create_parameter_df(lower_bound=-3,
                                     upper_bound=5)
    df['hierarchicalOptimization'] = 0
    df.parameterScale = 'log10'
    df.estimate = 1

    df.loc[df.index.str.startswith('scaling_'), 'hierarchicalOptimization'] = 1
    df.loc[df.index.str.startswith('sigma_'), 'hierarchicalOptimization'] = 1

    df.loc[df.index.str.startswith('offset_'), 'parameterScale'] = 'lin'
    df.loc[df.index.str.startswith('offset_'), 'lowerBound'] = -1e3
    df.loc[df.index.str.startswith('offset_'), 'upperBound'] = 1e3

    df.loc[df.index.str.startswith('offset_protein_total_genotypespecific_'),
           'hierarchicalOptimization'] = 0

    df.to_csv(parameter_file, sep="\t", index=True)


def update_ccle():
    condition_file_template = 'exp_table.tsv'
    measurement_file_template = 'measurements_Training_dataset1.tsv'
    sbml_file = 'CS_Signalling_ERBB_RAS_AKT_CCLE_petab.xml'
    condition_file = 'conditions_petab.tsv'
    measurement_file = 'measurements_CCLE_petab.tsv'
    parameter_file = 'parameters_CCLE_petab.tsv'

    sbml_reader = libsbml.SBMLReader()
    sbml_document = sbml_reader.readSBML(sbml_file)
    sbml_model = sbml_document.getModel()

    update_condition_table(condition_file_template, condition_file, sbml_model)
    update_measurement_table(measurement_file_template, measurement_file)
    create_parameter_table(sbml_file, condition_file, measurement_file,
                           parameter_file)

    # check for valid PEtab
    pp = petab.Problem.from_files(
        sbml_file=sbml_file,
        condition_file=condition_file,
        measurement_file=measurement_file,
        parameter_file=parameter_file)
    petab.lint_problem(pp)


def update_ccle_mclp():
    condition_file_template = 'exp_table.tsv'
    measurement_file_template = 'measurements_Training_dataset1-2.tsv'
    sbml_file = 'CS_Signalling_ERBB_RAS_AKT_CCLE_MCLP_petab.xml'
    condition_file = 'conditions_petab.tsv'
    measurement_file = 'measurements_CCLE_MCLP_petab.tsv'
    parameter_file = 'parameters_CCLE_MCLP_petab.tsv'

    sbml_reader = libsbml.SBMLReader()
    sbml_document = sbml_reader.readSBML(sbml_file)
    sbml_model = sbml_document.getModel()

    #update_condition_table(condition_file_template, condition_file, sbml_model)
    update_measurement_table(measurement_file_template, measurement_file)
    create_parameter_table(sbml_file, condition_file, measurement_file,
                           parameter_file)

    # check for valid PEtab
    pp = petab.Problem.from_files(
        sbml_file=sbml_file,
        condition_file=condition_file,
        measurement_file=measurement_file,
        parameter_file=parameter_file)
    petab.lint_problem(pp)

def main():
    # update_ccle()
    update_ccle_mclp()


if __name__ == '__main__':
    main()
