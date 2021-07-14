#!/usr/bin/env python
"""Initialize the project's data space.

Iterates over all defined state points and initializes
the associated job workspace directories.
The result of running this file is the creation of a signac workspace:
    - signac.rc file containing the project name
    - signac_statepoints.json summary for the entire workspace
    - workspace/ directory that contains a sub-directory of every individual statepoint
    - signac_statepoints.json within each individual statepoint sub-directory.
"""

import signac
import logging
from collections import OrderedDict
from itertools import product
import numpy as np

def get_parameters():
    '''
    Parameters:
    -----------

    System generation parameters:
    -----------------------------
    molecule : str
        Name of the molecule used to build the system.
        Must match one of the json files in uli-init/compounds
    para_weight : float; between 0 and 1
        The relative amount of para conformations in the system
        1 = All para, 0 = All meta
    density : float
        The density of the system in g/cm^3.
        PEEK and PEKK are both around 1.3 - 1.4 g/cm^3
    n_compounds : list
        A list of the number of molecules of a given length
        Must be the same legnth as polymer_lengths list(s)
        Corresponds to the number of specific molecules
        at the same index position in polymer_lengths
        See pdi parameter
    polymer_lengths : list
        A list of the number of monomer units in a single molecule
        Must be the same legnth as n_compounds list(s)
        See pdi parameter
    sample_pdi : bool
        Instruct uli-init to generate a distribution using a combination
        of pdi, Mn, Mw. This will override n_compound and polymer_length
        parameters
    pdi : float
        A PDI (poly-dispersity index) value of the generated system.
        pdi = Mn/Mw
    Mn : int
        The most frequent polymer length of a polydisperse system
        Used in conjunction with pdi to determine distribution
        of polymer lengths in the system
    Mw : int
        The weight average of the polymer distribution.
    forcefield : str options are 'gaff' or 'opls'
        The forcefield type to use when calling Foyer
    mass_dist : str
        Specify the distribution to be used when sampling from a pdi
        Options are: 'weibull' or 'gaussian'
    walls : bool
        If set to True, then hoomd will set LJ wall potentials on both
        sides of the volume along the x-axis.
        This should be True if you want to produce a trajectory with a flat
        surface.

    Simulation parameters:
    ----------------------


    ------------
    Other Notes:
    ------------
    All temperatures are entered as reduced temperature units

    If you want to sample from a PDI:
        Change the polymer length lines to [None]
        Change the n_compounds lines to [None]

    If you only want to run a quench simulation
        Comment out kT_anneal, anneal_sequence lines

    If you only want to run an anneal simulation
        Comment out kT_quench and n_steps lines

    Don't forget to change the name of the project
    project = signac.init_project("project-name")
    '''

    parameters = OrderedDict()
    ### System generation parameters ###
    parameters["system_type"] = [
	    "pack",
            #"stack",
            #"lamellar",
            #"coarse_grain",
            ]
    parameters["molecule"] = ['PEEK',
                             #'PEKK'
                             ]
    parameters["para_weight"] = [0.70]
    parameters["monomer_sequence"] = [None]
    parameters["density"] = [0.8]
    parameters["n_compounds"] = [[1]]
    parameters["polymer_lengths"] = [[50]]   
    parameters["pdi"] = [None]
    parameters["Mn"] = [None]
    parameters["Mw"] = [None]
    parameters['mass_dist'] = ['weibull']
    parameters["forcefield"] = ['gaff']
    parameters["remove_hydrogens"] = [True]
    parameters["system_seed"] = [24]
    parameters["kwarg"] = [{}]

    ### Simulation parameters ###
    parameters["tau_kt"] = [0.1]
    parameters["tau_p"] = [None]
    parameters["pressure"] = [None]
    parameters["dt"] = [0.001]
    parameters["e_factor"] = [0.5]
    parameters["sim_seed"] = [42]
    parameters["neighbor_list"] = ["cell"]
    parameters["walls"] = [False]
    parameters["shrink_kT"] = [0.2]
    parameters["shrink_steps"] = [5e6]
    parameters["procedure"] = [
            "quench",
            #"anneal"
        ]

    ### Quench related parameters ###
    parameters["kT_quench"] = [1.5]
    parameters["n_steps"] = [1e7]

    ### Anneal related parameters ###
    # List of [initial kT, final kT] Reduced Temps
    #parameters["kT_anneal"] = [
    #        [6.0, 2.0]
    #    ]     
    # List of lists of number of steps 
    #parameters["anneal_sequence"] = [
    #        [2e5, 1e5, 3e5, 5e5, 5e5, 1e5]
    #    ]
    #parameters["schedule"] = [None]

    ### Do a few checks to save headaches later:
    if [None] not in [
			parameters["para_weight"], parameters["monomer_sequence"]
			]:
        raise ValueError(
                "You can only use one of `para_weight` and "
                "`monomer_sequence`. Set the other to [None]"
                )
    return list(parameters.keys()), list(product(*parameters.values()))

custom_job_doc = {} # added keys and values to be added to each job document created
                    # leave blank to create for job doc entries
def main():
    project = signac.init_project("project")
    param_names, param_combinations = get_parameters()
    # Create the generate jobs
    for params in param_combinations:
        parent_statepoint = dict(zip(param_names, params))
        parent_job = project.open_job(parent_statepoint)
        parent_job.init()
        try:
            parent_job.doc.setdefault("steps", parent_statepoint["n_steps"])
        except:
            parent_job.doc.setdefault("steps", np.sum(parent_statepoint["anneal_sequence"]))
            parent_job.doc.setdefault("step_sequence", parent_statepoint["anneal_sequence"])
        if any([parent_job.sp['Mn'], parent_job.sp['pdi'], parent_job.sp['Mw']]):
            parent_job.doc.setdefault("sample_pdi", True)
        else:
            parent_job.doc.setdefault("sample_pdi", False)
    
    if custom_job_doc:
        for key in custom_job_doc:
            parent_job.doc.setdefault(key, custom_job_doc[key])

    project.write_statepoints()

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
