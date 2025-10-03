import sys                                  #Have to code for multiple seeds,have to add custom cylindrical camkii_rel and just run for diff ca,diff configuration of camkii??
import os
import random
import copy
import numpy as np


# set recursively region to all nodes of type m.ExprNodeType.LEAF
def set_region_to_all_count_term_leaf_nodes(count_term_node, region):
    if not count_term_node:
        return

    if count_term_node.node_type != m.ExprNodeType.LEAF:
        set_region_to_all_count_term_leaf_nodes(count_term_node.left_node, region)
        set_region_to_all_count_term_leaf_nodes(count_term_node.right_node, region)
    else:
        count_term_node.region = region


def load_bngl_observables_and_create_psd_and_spine_variants(model, file_name, spine, psd, sampling_periodicity, seed):

    # create a helper object just to load the observables from BNGL,
    observables = m.Observables()
    observables.load_bngl_observables(file_name)

    for base_count in observables.counts:
        # count in psd
        count_psd = copy.deepcopy(base_count)
        count_psd.name = str(count_psd.name + '_psd') 
        count_psd.file_name = str('./react_data/seed_' + str(seed).zfill(5) + '/' + count_psd.name + '.dat')
        count_psd.every_n_timesteps = sampling_periodicity
        set_region_to_all_count_term_leaf_nodes(count_psd.expression, psd)
        model.add_count(count_psd)

        # count everything that is in spine but not in psd
        count_reg_spine_minus_psd = copy.deepcopy(base_count)
        count_reg_spine_minus_psd.name = str(count_reg_spine_minus_psd.name + '_reg_spine_minus_psd') 
        count_reg_spine_minus_psd.file_name = str('./react_data/seed_' + str(seed).zfill(5) + '/' + count_reg_spine_minus_psd.name + '.dat')
        count_reg_spine_minus_psd.every_n_timesteps = sampling_periodicity
        set_region_to_all_count_term_leaf_nodes(count_reg_spine_minus_psd.expression, spine - psd)
        model.add_count(count_reg_spine_minus_psd)


MCELL_PATH=os.environ.get('MCELL_PATH','')
print(MCELL_PATH)
if MCELL_PATH:
    lib_path=os.path.join(MCELL_PATH,'lib')
    sys.path.append(lib_path)
    print('Yay!')
else:
    print("Error")
    sys.exit(1)

import mcell as m 
if __name__ == "__main__":


    SEED = 2
    n_pulses = 1 #get_parameter_value(n_pulses)
    dt = 0.0001
    
    V_true = 0.018822e-15
    NA = 6.022e23/1e6 #1e-6


    tauR = 0.002 #time constant for Ca decay
    tauF = 0.01

    t1 = 0
    t0 = 0.01
    x = np.array(np.arange(0,0.08,dt))
    alpha = 2000*(NA*V_true)*((x-t1)/tauR)*np.exp(-(x-t1)/tauF)*dt
    params = m.bngl_utils.load_bngl_parameters('ryr2_CaM.bngl')
    TIME_STEP = 1e-6
    ITERATIONS = int(params['ITERATIONS'])
                      
    species_dict = dict([
    ("Ca", NA*V_true*0.1),
    ("CaM(C~0,N~0,ryr2)", NA*V_true*4),#NA*V_true*30,
    ("RyR2(r!1,l!4,Ca~0,cam,S~0).RyR2(r!2,l!1,Ca~0,cam,S~0).RyR2(r!3,l!2,Ca~0,cam,S~0).RyR2(r!4,l!3,Ca~0,cam,S~0)", 4)])
    seed = 10
    
    import spine133_geometry as geometry
    model = m.Model()


    spine = geometry.d000p_sp133_cut
    model.add_geometry_object(spine)
    PSD = geometry.d000p_sp133_PSD
    model.add_geometry_object(PSD)
   
    



    for name in species_dict.keys():
        print(name)
        rel_site = m.ReleaseSite(
            name = 'Release_of_' + name,
            complex = m.Complex(name),
            region = spine,
            number_to_release = (int)(species_dict[name])
        )
        model.add_release_site(rel_site)
    model.load_bngl("ryr2_CaM.bngl",'./react_data_wm/seed_' + str(SEED).zfill(5) + '/', spine)
    for c in model.counts:
        c.every_n_timesteps = ITERATIONS/100000
    # for rel in range(len(x)):
    #     if alpha[rel] > 0:
    #         new_rel = m.ReleaseSite(
    #             name = 'rel_%i' % rel,
    #             complex = m.Complex('Ca'),
    #             ## This location for spine048
    #             #location = [0.009,-0.045,-0.1],
    #             ## This location for spine133
    #             location = [5.868,5.62319,5.88514],
    #             number_to_release =alpha[rel],
    #             release_time = t0+x[rel]
    #     )
    #         model.add_release_site(new_rel)

    tp = m.SurfaceClass(
        name = 'transparent',
        type = m.SurfacePropertyType.TRANSPARENT,
        affected_complex_pattern = m.AllMolecules
    )

    model.add_surface_class(tp)
    PSD.surface_class = tp
    

    model.config.time_step = TIME_STEP
    model.config.seed = SEED
    model.config.total_iterations = ITERATIONS
    model.config.partition_dimension = 20#1
    model.config.subpartition_dimension = 0.05

    load_bngl_observables_and_create_psd_and_spine_variants(
        model, 'ryr2_CaM.bngl', spine, PSD, ITERATIONS/100000, SEED) 

    viz_output = model.VizOutput(
        mode = m.VizMode.ASCII,
        output_files_prefix = './viz_data/seed_' + str(SEED).zfill(5) + '/Scene',
        every_n_timesteps = SAMPLING_PERIODICITY
    )

    model.initialize()

    model.export_data_model()
    model.run_iterations(ITERATIONS)
    model.end_simulation()
