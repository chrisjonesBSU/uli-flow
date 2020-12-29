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
    print('in sampled function')
    print(current_step(job))
    print(job.doc.steps)
    already_sampled = current_step(job) >= job.doc.steps
    print(already_sampled)
    return current_step(job) >= job.doc.steps


@MyProject.label
def initialized(job):
    return job.isfile("init.pdb")


@directives(executable="python -u")
@directives(ngpu=1)
@MyProject.operation
# @MyProject.pre.after(initialize)
@MyProject.post(sampled)
def sample(job):
    print('SAMPLED FUNCTION OUTPUT:')
    sampled_output = sampled(job)
    print(sampled_output)
    from uli_init import simulate
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
                pdi = job.sp['pdi'],
                M_n = job.sp['M_n'],
                remove_hydrogens = job.sp['remove_hydrogens'],
                seed = job.sp['system_seed']
            )

        #system.system_mb.save('init.pdb')
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
                gsd_write = 1e5,
                log_write = 1e4
                )
        logging.info("Simulation object generated...")

        if job.sp['procedure'] == "quench":
            logging.info("Beginning quench simulation...")
            simulation.quench(
                    kT = job.sp['kT_quench'],
                    n_steps = job.sp['n_steps'],
                    shrink_kT = 10,
                    shrink_steps = 1e6
                    )

        elif job.sp['procedure'] == "anneal":
            logging.info("Beginning anneal simulation...")
            simulation.anneal(
                    kT_init = job.sp['kT_anneal'][0],
                    kT_final = job.sp['kT_anneal'][1],
                    step_sequence = job.sp['anneal_sequence'],
                    schedule = job.sp['schedule'],
                    shrink_kT = 10,
                    shrink_steps = 1e6
                    )


if __name__ == "__main__":
    MyProject().main()
