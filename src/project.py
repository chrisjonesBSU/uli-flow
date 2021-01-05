"""Define the project's workflow logic and operation functions.

Execute this script directly from the command line, to view your project's
status, execute operations and submit them to a cluster. See also:

    $ python src/project.py --help
"""
from flow import FlowProject, directives
from flow.environment import DefaultSlurmEnvironment
from flow.environments.xsede import BridgesEnvironment, CometEnvironment

class MyProject(FlowProject):
    pass

class Fry(DefaultSlurmEnvironment):
    hostname_pattern = "fry"
    template = "fry.sh"

    @classmethod
    def add_args(cls, parser):
        parser.add_argument(
            "--partition", default="batch", help="Specify the partition to submit to."
        )


# Definition of project-related labels (classification)
def current_step(job):
    import gsd.hoomd

    if job.isfile("sim_traj.gsd"):
        with gsd.hoomd.open(job.fn("sim_traj.gsd")) as traj:
            return traj[-1].configuration.step
    return -1


@MyProject.label
def sampled(job):
    return current_step(job) >= job.doc.steps


@MyProject.label
def initialized(job):
    return job.isfile("init.pdb")

@MyProject.label
def rdf_done(job):
    return job.isfile("rdf.csv")

@MyProject.label
def ind_sampling_done(job):
    return job.isfile("sim_traj_equil.log")

@directives(executable="python -u")
@directives(ngpu=1)
@MyProject.operation
# @MyProject.pre.after(initialize)
@MyProject.post(sampled)
def sample(job):
    from uli_init import simulate
    from uli_init.utils import base_units, unit_conversions
    import numpy as np
    import os
    import logging

    with job:
        logging.info("Creating system...")
        system = simulate.System(
                molecule = job.sp['molecule'],
                para_weight = job.sp['para_weight'],
                density = job.sp['density'],
                n_compounds = job.sp['n_compounds'],
                polymer_lengths = job.sp['polymer_lengths'],
                forcefield = job.sp['forcefield'],
                sample_pdi = job.doc.sample_pdi,
                pdi = job.sp['pdi'],
                Mn = job.sp['Mn'],
                Mw = job.sp['Mw'],
                mass_dist_type = job.sp['mass_dist'],
                remove_hydrogens = job.sp['remove_hydrogens'],
                seed = job.sp['system_seed']
            )

        if system.system_pmd:
            system.system_pmd.save('init.pdb', overwrite=True)
        else:
            system.system_mb.save('init.pdb', overwrite=True)

        job.doc['num_para'] = system.para
        job.doc['num_meta'] = system.meta
        job.doc['num_compounds'] = system.n_compounds
        job.doc['polymer_lengths'] = system.polymer_lengths

        logging.info("System generated...")
        logging.info("Starting simulation...")

        simulation = simulate.Simulation(
                system,
                target_box = None,
                r_cut = 1.2,
                e_factor = job.sp['e_factor'],
                tau = job.sp['tau'],
                dt = job.sp['dt'],
                seed = job.sp['sim_seed'],
                auto_scale = True,
                ref_units = None,
                mode = "gpu",
                gsd_write = max([int(job.doc['steps']/100), 1]),
                log_write = max([int(job.doc['steps']/10000), 1])
                )

        logging.info("Simulation object generated...")
        job.doc['ref_energy'] = simulation.ref_energy
        job.doc['ref_distance'] = simulation.ref_distance
        job.doc['ref_mass'] = simulation.ref_mass
        job.doc['real_timestep'] = unit_conversions.convert_to_real_time(simulation.dt,
                                                    simulation.ref_energy,
                                                    simulation.ref_distance,
                                                    simulation.ref_mass)
        job.doc['time_unit'] = 'fs'
        job.doc['steps_per_frame'] = simulation.gsd_write
        job.doc['steps_per_log'] = simulation.log_write

        if job.sp['procedure'] == "quench":
            job.doc['T_SI'] = unit_conversions.kelvin_from_reduced(job.sp['kT_quench'],
                                                    simulation.ref_energy)
            job.doc['T_unit'] = 'K'
            logging.info("Beginning quench simulation...")
            simulation.quench(
                    kT = job.sp['kT_quench'],
                    n_steps = job.sp['n_steps'],
                    shrink_kT = 10,
                    shrink_steps = 1e6
                    )

        elif job.sp['procedure'] == "anneal":
            logging.info("Beginning anneal simulation...")
            if not job.sp['schedule']:
                kT_list = np.linspace(job.sp['kT_anneal'][0],
                                      job.sp['kT_anneal'][1],
                                      len(job.sp['anneal_sequence']),
                                      )
                kT_SI = [unit_conversions.kelvin_from_reduced(kT, simulation.ref_energy)
                            for kT in kT_list]
                job.doc['T_SI'] = kT_SI
                job.doc['T_unit'] = 'K'

            simulation.anneal(
                    kT_init = job.sp['kT_anneal'][0],
                    kT_final = job.sp['kT_anneal'][1],
                    step_sequence = job.sp['anneal_sequence'],
                    schedule = job.sp['schedule'],
                    shrink_kT = 10,
                    shrink_steps = 1e6
                    )

@directives()
@MyProject.operation
#@MyProject.pre.after(sampled)
@MyProject.post(rdf_done)
#@MyProject.with_job
def post_process(job):
    '''
    X 1. Independence sampling using .log file
        - Update job doc
        - Save the sampled data to a new .log file
    2. Compute average RDF over 10-15 frames
        - Check a few different atom types
        - Save results to .csv files
    3. Compute MSD
        - Check a few different atom types
        - Save results to .csv files
    4. Save some plots (PE, RDF, MSD) in pdf format
    5. 

    '''
    import gsd.hoomd
    import freud
    import numpy as np
    from cme_lab_utils import msd, rdf, sampler, gsd_utils
    import os
    import logging
    import matplotlib.pyplot as plt
    
    # Perform independence sampling:
    if job.sp['procedure'] == 'quench':
        start_index = 0
    elif job.sp['procedure'] == 'anneal': # Only want to sample from last temp
            start_index = -job.sp['anneal_sequence'][-1] // job.doc['steps_per_log']

    job_log_file = np.genfromtxt(job.fn('sim_traj.log'),
                                delimiter='\t',
                                names=True
                                )
    pe = job_log_file['potential_energy']
    samples = sampler.Mbar(pe[start_index:], nskip=1)
    job.doc['production_start'] = samples.start # Starting index of production region
    job.doc['production_ineff'] = samples.ineff # Statistical inefficiency of prod region
    job.doc['production_size'] = samples.production_size # Size of prod region
    job.doc['num_ind_samples'] = len(samples.indices)
    sampled_data = job_log_file[samples.start:][samples.indices]
    col_names = [name for name in job_log_file.dtype.names]
    headers = "{}\t"*(len(col_names) - 1)+"{}"
    np.savetxt(job.fn('sim_traj_equil.log'),
                sampled_data,
                header = headers.format(*col_names)
              )
    logging.info("Finished independence sampling...") 

    # Calculate some RDFs and MSDs from GSD files and save results to txt files
    types = [['ca', 'ca'], ['oh', 'oh']]
    rdf_dir = os.path.join(job.ws, 'rdf-results')
    if not os.path.exists(rdf_dir):
        os.mkdir(rdf_dir)
    for pair in types:
        pair_rdf = rdf.gsd_rdf(job.fn("sim_traj.gsd"),
                               pair[0],
                               pair[1],
                               start=-10
                               )
        x = pair_rdf.bin_centers
        y = pair_rdf.rdf
        fig = plt.figure()
        plt.plot(x, y)
        plt.title("RDF: Atom Types {} {}".format(*pair))
        plt.xlabel('r')
        plt.ylabel('g(r)')
        fig.savefig(os.path.join(rdf_dir, "{}_{}.pdf".format(*pair)))
        np.savetxt(os.path.join(rdf_dir, '{}_{}.txt'.format(*pair)),
                   np.transpose([x, y]),
                   header = "x,y", 
                   delimiter=","
                  )
    logging.info("Finished RDF calculations...")

    # Start MSD calculations
    types = ["ca", "oh"]
    msd_dir = os.path.join(job.ws, 'msd-results')
    if not os.path.exists(msd_dir):
        os.mkdir(msd_dir)
    for atom_type in types:
        msd_results = msd.msd_from_gsd(job.fn("sim_traj.gsd"),
                                       start=-10, stop=-1,
                                       atom_type = atom_type,
                                       msd_mode="window"
                                      )
        
        seconds = np.arange(0, len(msd_results), 1) * job.doc['real_timestep'] * job.doc['steps_per_frame']
        fig = plt.figure()
        plt.plot(seconds, msd_results)
        plt.title("Mean Square Displacement")
        plt.xlabel("Seconds")
        plt.ylabel("MSD ($\AA^2$)")
        fig.savefig(os.path.join(msd_dir, "{}.pdf".format(atom_type)))
        np.savetxt(os.path.join(msd_dir, '{}.txt'.format(atom_type)),
                   msd_results
                   )
    logging.info("Finished MSD Calculations")


                          







    


                           



if __name__ == "__main__":
    MyProject().main()
