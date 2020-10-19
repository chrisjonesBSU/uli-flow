"""Define the project's workflow logic and operation functions.

Execute this script directly from the command line, to view your project's
status, execute operations and submit them to a cluster. See also:

    $ python src/project.py --help
"""
from flow import FlowProject, directives
from flow.environment import DefaultSlurmEnvironment
from flow.environments.xsede import BridgesEnvironment, CometEnvironment


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

    if job.isfile("trajectory.gsd"):
        with gsd.hoomd.open(job.fn("trajectory.gsd")) as traj:
            return traj[-1].configuration.step
    return -1


@MyProject.label
def sampled(job):
    return current_step(job) >= job.doc.steps


@MyProject.label
def initialized(job):
    return job.isfile("init.hoomdxml")


@directives(executable="python -u")
@directives(ngpu=1)
@MyProject.operation
# @MyProject.pre.after(initialize)
@MyProject.post(sampled)
def sample(job):
    from uli_init import simulate
    import os
    import logging

    with job:
        
        system = simulate.System(
                molecule = job.sp['molecule'],
                para_weight = job.sp['para_weight'],
                density = job.sp['density'],
                n_compounds = job.sp['n_compounds'],
                polymer_lengths = job.sp['polymer_lengths'],
                forcefield = job.sp['forcefield'],
                pdi = job.sp['pdi'],
                M_n = job.sp['M_n'],
                remove_hydrogens = job.sp['remove_hydrogens']
            )

        simulation = simulate.Simulation(
                system,
                target_box = None,
                r_cut = 1.2,
                e_factor = job.sp['e_factor'],
                tau = job.sp['tau'],
                dt = job.sp['dt'],
                auto_scale = True,
                ref_units = None,
                mode = "gpu",
                gsd_write = 1e5,
                log_write = 1e3
                )

        if job.sp['procedure'] == "quench":
            simulation.quench(
                    kT = job.sp['kT_quench'],
                    n_steps = job.sp['n_steps'],
                    shrink_kT = 10,
                    shrink_n_steps = 1e6
                    )

        elif job.sp['procedure'] == "anneal":
            simulation.anneal(
                    kT_init = job.sp['kT_anneal'][0],
                    kT_final = job.sp['kT_anneal'][1],
                    step_sequence = job.sp['step_sequence'],
                    schedule = job.sp['schedule'],
                    shrink_kT = 10,
                    shrink_n_steps = 1e6
                    )


