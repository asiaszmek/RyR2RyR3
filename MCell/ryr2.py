import sys                                  #Have to code for multiple seeds,have to add custom cylindrical camkii_rel and just run for diff ca,diff configuration of camkii??
import os
import random
MCELL_PATH=os.environ.get('MCELL_PATH','')
if MCELL_PATH:
    lib_path=os.path.join(MCELL_PATH,'lib')
    sys.path.append(lib_path)
   # print('yessuh')
else:
    print("Error")
    sys.exit(1)

import mcell as m 
for x in range(20):
        # viz_output = m.VizOutput(
        #     output_files_prefix = './camkii_final1/viz_data/seed_00001/Scene',
# )
        # random.seed(x)
    box_vertex_list = [
            [ 0.025000000372529, -0.25, -0.25 ],
            [ 0.025000000372529, 0.25, -0.25 ],
            [ 0.524999976158142, 0.25, 0.25 ],
            [ 0.524999976158142, 0.25, -0.25 ],
            [ 0.025000000372529, 0.25, 0.25 ],
            [ 0.524999976158142, -0.25, 0.25 ],
            [ 0.025000000372529, -0.25, 0.25 ],
            [ 0.524999976158142, -0.25, -0.25 ]
        
        ] 

    box_wall_list = [
            # [1, 2, 0] defines a triangle connecting vertices 
            # [-0.1, -0.1, 0.1], [-0.1, 0.1, -0.1], and
            # [-0.1, -0.1, -0.1]
            [ 7, 0, 3 ],
            [ 3, 0, 1 ],
            [ 4, 1, 6 ],
            [ 0, 6, 1 ],
            [ 1, 4, 3 ],
            [ 2, 3, 4 ],
            [ 7, 3, 5 ],
            [ 2, 5, 3 ],
            [ 6, 5, 4 ],
            [ 2, 4, 5 ],
            [ 5, 0, 7 ],
            [ 6, 0, 5 ]
        ] 
    postsynaptic_PSD_wall_indices = [
            2, 3
        ]
    postsynaptic_PSD = m.SurfaceRegion(
            name = 'PSD',
            wall_indices = postsynaptic_PSD_wall_indices
        )

    box_vertex_list2 = [
            [ 0.0875197052955627, -0.150000005960464, -0.150000005960464 ],
            [ 0.0875197052955627, -0.150000005960464, 0.150000005960464 ],
            [ 0.0875197052955627, 0.150000005960464, -0.150000005960464 ],
            [ 0.0875197052955627, 0.150000005960464, 0.150000005960464 ],
            [ 0.387519717216492, -0.150000005960464, -0.150000005960464 ],
            [ 0.387519717216492, -0.150000005960464, 0.150000005960464 ],
            [ 0.387519717216492, 0.150000005960464, -0.150000005960464 ],
            [ 0.387519717216492, 0.150000005960464, 0.150000005960464 ]
        
        ] 

    box_wall_list2 = [
            # [1, 2, 0] defines a triangle connecting vertices 
            # [-0.1, -0.1, 0.1], [-0.1, 0.1, -0.1], and
            # [-0.1, -0.1, -0.1]
            [ 0, 1, 3 ],
            [ 2, 3, 7 ],
            [ 6, 7, 5 ],
            [ 4, 5, 1 ],
            [ 2, 6, 4 ],
            [ 7, 3, 1 ]
            
        ] 
    

    box = m.GeometryObject(
            name = 'box',
            vertex_list = box_vertex_list,
            wall_list = box_wall_list,
            surface_regions = [postsynaptic_PSD]
        )
    box2 =m.GeometryObject(
            name = 'box2',
            vertex_list = box_vertex_list2,
            wall_list = box_wall_list2,
        ) 
    Ca_pat = m.ReleasePattern(
            name = 'Ca_pat',
            release_interval = 100e-3,
            train_duration = 160e-3,
            train_interval=0.5,
            number_of_trains = 1
        )
    Ca_waveform = m.ReactionRule(
    name = 'Ca_waveform',
    reactants = [ m.Complex('Ca') ],
    products = [ ],
    fwd_rate = 50
    )
    trans_a = m.SurfaceClass(
    name = 'trans_a',
    type = m.SurfacePropertyType.TRANSPARENT,
    affected_complex_pattern = m.AllMolecules
    )
    # box2.surface_class = trans_a
    # box.is_bngl_compartment = True
    # box.surface_compartment_name='O1M'
    release_site_a = m.ReleaseSite(
            name = 'rel_a', 
            complex = m.Complex('CaM(C~0,N~0,CaMKII)'),# compartment_name='box'),
            region=box2,
            #location=(0,0,0),
            number_to_release = 600      #change
        )
    release_site_b = m.ReleaseSite(
            name = 'rel_b', 
            complex = m.Complex('Ca()'),#compartment_name='box'), 
            release_pattern=Ca_pat,
            region=box, 
            number_to_release =4320         #change
        )
    release_site_c = m.ReleaseSite(
            name = 'rel_c', 
            complex = m.Complex('CaMKII(r!1,l!6,T286~U,CaM).CaMKII(r!2,l!1,T286~U,CaM).CaMKII(r!3,l!2,T286~U,CaM).CaMKII(r!4,l!3,T286~U,CaM).CaMKII(r!5,l!4,T286~U,CaM).CaMKII(r!6,l!5,T286~U,CaM)'),# compartment_name='box'),
            region=box2,
            #location=(0,0,0),
            number_to_release = 172        #change
        )
    
    SEED=random.randint(1,10000)
    model = m.Model()
    model.config.seed = SEED
    model.add_reaction_rule(Ca_waveform)
        # model.add_viz_output(viz_output)  
    model.add_geometry_object(box)
    model.add_geometry_object(box2)
    model.add_surface_class(trans_a)
    model.add_release_site(release_site_a)
    model.add_release_site(release_site_b)
    model.add_release_site(release_site_c)



    MODEL_PATH = os.path.dirname(os.path.abspath(__file__))


    model.load_bngl(
            file_name = os.path.join(MODEL_PATH, 'ryr2_CaM.bngl'),
            observables_path_or_file = './camkii/100ms6/seed_0000{}/'.format(x)
        )    
    ITERATIONS = 150000
    model.config.total_iterations = ITERATIONS 

    model.initialize()
        # model.export_viz_data_model()
    model.run_iterations(ITERATIONS)
    model.end_simulation()










