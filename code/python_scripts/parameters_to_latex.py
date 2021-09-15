## Imports

import pandas as pd


def create_parameter_table(datapath):
    # Get and clean up parameters
    parameters = pd.read_csv(os.path.join(datapath, 'sim_shared_parameters.csv'), index_col=False)
    parameters = parameters.drop(columns = ['layer_descr', 'model_names_1', 'model_names_2', 'model_names_3'])

    # Clean up

    parameters = parameters.rename(columns={"taskweights_1":"taskweights_aeh",
                               "taskweights_2":"taskweights_kon",
                               "taskweights_3":"taskweights_sum",
                               "phaseweights_1": "phaseweights_cue",
                               "phaseweights_2": "phaseweights_trial",
                               "phaseweights_3":"phaseweights_rest",
                               "response_1": "response_cue",
                               "response_2": "response_trial",
                               "response_3":"response_rest"
                               })


    # Edit parameter table

    # Save parameters as LateX table
    datapath = '/Users/leonardo.pettini/Desktop/Masterarbeit_BCCN/ROI_task_positive_network/group'
    parametersLateX = parameters.to_latex(index=False)
    tableFileName = 'sim_shared_parameters.tex'
    parametersLateXPath =  os.path.join(datapath, tableFileName)

    with open(parametersLateXPath,'w') as tf:
        tf.write(parametersLateX)

    print(parametersLateX)

        return parametersLateX