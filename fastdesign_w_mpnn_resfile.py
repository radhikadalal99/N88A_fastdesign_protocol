#!/usr/bin/python
from pyrosetta import *
import sys
init()

itr=sys.argv[1]

xml_string='''
 <SCOREFXNS>
    <ScoreFunction name="ref2015" weights="ref2015"/>
    <ScoreFunction name="sfxn_design" weights="ref2015" >
        <Reweight scoretype="approximate_buried_unsat_penalty" weight="5" />
        <Set approximate_buried_unsat_penalty_burial_probe_radius="2.3" />
        <Set approximate_buried_unsat_penalty_burial_atomic_depth="4.0" />
        <Set approximate_buried_unsat_penalty_hbond_energy_threshold="-0.5" />
    </ScoreFunction>
  </SCOREFXNS>
'''
xml = pyrosetta.rosetta.protocols.rosetta_scripts.XmlObjects.create_from_string(xml_string)
sfxn_design=xml.get_score_function("sfxn_design")


score_manager = pyrosetta.rosetta.core.scoring.ScoreTypeManager()
score_term = score_manager.score_type_from_name("hbond_sc")
sfxn_design.set_weight(score_term, 2.0)
score_term1 = score_manager.score_type_from_name("res_type_constraint")
sfxn_design.set_weight(score_term1,1.5)

p=pose_from_file("relaxed_recleaned_7kx0_fab.pdb")

#pyrosetta.rosetta.protocols.protein_interface_design.FavorNativeResidue(p,1.5)
restypecst=pyrosetta.rosetta.protocols.constraint_generator.ResidueTypeConstraintGenerator()
restypecst.set_favor_native_bonus(1.5)
restypecst.apply(p)


tf = pyrosetta.rosetta.core.pack.task.TaskFactory()
tf.push_back(pyrosetta.rosetta.core.pack.task.operation.IncludeCurrent())
tf.push_back(pyrosetta.rosetta.core.pack.task.operation.NoRepackDisulfides())
# Include the resfile
tf.push_back(pyrosetta.rosetta.core.pack.task.operation.ReadResfile("./mpnn_based_resfile.txt"))
# Convert the task factory into a PackerTask
packer_task=tf.create_task_and_apply_taskoperations(p)
fast_design = pyrosetta.rosetta.protocols.denovo_design.movers.FastDesign(scorefxn_in=sfxn_design, standard_repeats=5)
fast_design.cartesian(False)
fast_design.set_task_factory(tf)
#mm=pyrosetta.rosetta.core.kinematics.MoveMap()
#mm.init_from_file("./movemap.txt")
#fast_design.set_movemap_factory(movemap)
fast_design.min_type("lbfgs_armijo_nonmonotone")

fast_design.apply(p)


scorefxn=pyrosetta.create_score_function("ref2015")
relax=pyrosetta.rosetta.protocols.relax.FastRelax()
relax.set_scorefxn(scorefxn)
relax.constrain_relax_to_start_coords(True) #maintain bb positions
relax.ramp_down_constraints(True)
relax.constrain_coords(True)
relax.set_task_factory(tf)
for i in range(1,6):
        relax.apply(p)

p.dump_pdb("%s_fd_redesign.pdb"%itr)
