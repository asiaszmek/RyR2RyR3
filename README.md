Stochastic model of ryanodine type 2 and 3 receptor activation by calcium and gating by calmodulin
=======================


Introduction
-----------------------
This repo contains two kinetic models: one of the ryanodine type 2 receptor (RyR2) inactivation by calmodulin and activation by calcium and ryanodine type 3 receptor (RyR3) activation by calcium and calmodulin. RyR2 is expressed mostly in cardiac myocites and has been quite thoroughly investigated. However, there is much less known about the RyR3, which is express throught the body. We used the known measurements of affinity of calmodulin binding to CaMBD2 domain of RyR2 to built both RyR2-CaM and RyR3-CaM [1, 2, 3] and dissociation rates [4, 5], because CaMBD2 is highly conserved for both ryanodine receptor types [6].

We chose a 4 step kinetic gating scheme of the Hodgkin-Huxley potassium channel to model calcium bindings to RyR (one calcium ion binding each subunit of the RyR tetramer) with Ca4-RyR switching to 2 open states (O1 and O2). Both O1 and O2 can switch to an additional closed state Ca4-RyR-C, which can also transition to an inactive state Ca4-RyR-I, which is a kinetic gating scheme proposed by [3].

Models have been fitted to experimental data. Calcium activation of RyR2 unbound to calmodulin was fitted to RyR2 opening time course from [7], and verfied by evaluating the fit to open probability, mean open and closed time provided by [3]. Calcium activation of CaM-bound RyR2 open probability was fitted to appropriate data from [3] and the mean open and closed time was validated to best match data from [3]. 

RyR3-CaM model was contructed based on the RyR3 model. Calcium binding rates to RyR3 and opening and closing rates of RyR3 (unbound to CaM) were fitted to data from [8]. The opening and closing rates of CaM-bound RyR3 were fitted to open probability in the function of calcium concentration combined from [8] and [9].

For RyR and RyR-CaM models we assumed that transition rates between Ca4-RyR and Ca4-RyR-O2 are at least two orders of magnitude lower than Ca4-RyR and Ca4-RyR-O1.

Models in the model directory are implemented using NeuroRD, a stochastic reaction-diffusion simulator implementing Gillespie algorith with asynchronous tau-leaping [10] (https://github.com/neurord/stochdiff). Models in the MCell directory are implemented in the MCell (https://mcell.org) environment using bngl (BioNetGen Language, https://bionetgen.org/). Kinetic constants of the NeuroRD model were were fitted using ajustador [11] (https://github.com/neurord/ajustador) and neurord_fit (https://github.com/neurord/neurord_fit).

Figures from journal articles investigating RyR2 and RyR3 properties and data files used to fit RyR2 and RyR3 models can be found in datasets_for_fitting directory. Models were fitted using following scripts: neurord_fit_Ca_RyR2.py,neurord_fit_Ca_RyR2CaM_steady_state.py,  neurord_fit_Ca_RyR2_steady_state.py, neurord_fit_Ca_RyR3_steady_state.py, neurord_fit_Ca_RyR3CaM_steady_state.py. For each part of the receptor gating kinetic scheme validation was performed by reproducing expreriments investigating probability, mean open and closed time in the function of calcium concentration. Scripts for validating fits can be found in following directories: fit_RyR2_Ca, fit_RyR2CaM_Ca, fit_RyR3_Ca, fit_RyR3CaM_Ca.

Models and fits were developed using git version control system with commits serving as log entries for the lab notebook.

[1] Lau K, Chan MMY, Van Petegem F. Lobe-Speciﬁc Calmodulin Binding to Different Ryanodine Receptor Isoforms. Biochemistry. 2014; 53(5):932–946. https://doi.org/10.1021/bi401502x, doi: 10.1021/bi401502x, pMID:
24447242.

[2] Søndergaard MT, Liu Y, Guo W, Wei J, Wang R, Brohus M, Overgaard MT, Chen SRW. Role of cardiac ryanodine receptor calmodulin-binding domains in mediating the action of arrhythmogenic calmodulin N-domain mutation N54I. The FEBS Journal. 2020; 287(11):2256–2280. https://febs.onlinelibrary.wiley.com/doi/abs/10.1111/febs.15147, doi: https://doi.org/10.1111/febs.15147.

[3] Xu L, Meissner G. Mechanism of Calmodulin Inhibition of Cardiac Sarcoplasmic Reticulum Ca2+ Release Channel (Ryanodine Receptor). Biophysical Journal. 2004; 86(2):797–804. https://www.sciencedirect.com/science/article/pii/S0006349504741557, doi: https://doi.org/10.1016/S0006-3495(04)74155-7.

[4] Balshaw DM, Xu L, Yamaguchi N, Pasek DA, Meissner G. Calmodulin Binding and Inhibition of Cardiac Muscle Calcium Release Channel (Ryanodine Receptor)*.
Journal of Biological Chemistry. 2001; 276(23):20144–20153. https://www.sciencedirect.com/science/article/pii/S0021925819404547, doi:https://doi.org/10.1074/jbc.M010771200.

[5] Wu X, Bers DM. Free and bound intracellular calmodulin measurements in cardiac myocytes. Cell Calcium. 2007; 41(4):353–364. https://www.sciencedirect.com/science/article/pii/S0143416006001618, doi:https://doi.org/10.1016/j.ceca.2006.07.011.

[6] Brohus M, Søndergaard MT, Wayne Chen SR, van Petegem F, Overgaard MT. Ca2+-dependent calmodulin binding to cardiac ryanodine receptor (RyR2) calmodulin-binding domains. Biochemical Journal. 2019 01; 476(2):193–209. https://doi.org/10.1042/BCJ20180545, doi: 10.1042/BCJ20180545.

[7] Zahradníková A, Zahradník I, Györke I, Györke S. Rapid Activation of the Cardiac Ryanodine Receptor by Submillisecond Calcium Stimuli. Journal of General Physiology. 1999 11; 114(6):787–798. https://doi.org/10.1085/
jgp.114.6.787, doi: 10.1085/jgp.114.6.787.

[8] Chen SRW, Li X, Ebisawa K, Zhang L. Functional Characterization of the Recombinant Type 3 Ca2+ Release Channel (Ryanodine Receptor) Expressed in HEK293 Cells*. Journal of Biological Chemistry. 1997; 272(39):24234–24246. https://www.sciencedirect.com/science/article/pii/S0021925819635891, doi:https://doi.org/10.1074/jbc.272.39.24234.

[9] Yamaguchi N, Xu L, Pasek DA, Evans KE, Chen SRW, Meissner G. Calmodulin Regulation and Identiﬁcation of Calmodulin Binding Region of Type-3 Ryanodine Receptor Calcium Release Channel. Biochemistry. 2005; 44(45):15074–15081. https://doi.org/10.1021/bi051251t, doi: 10.1021/bi051251t, pMID: 16274254.

[10] Jędrzejewski-Szmek Z, Blackwell KT. Asynchronous 𝜏-leaping. The Journal of Chemical Physics. 2016 03; 144(12):125104. https://doi.org/10.1063/1.4944575, doi: 10.1063/1.4944575.

[11] Jędrzejewski-Szmek Z, Abrahao KP, Jędrzejewska-Szmek J, Lovinger DM, Blackwell KT. Parameter Optimization Using Covariance Matrix Adaptation—Evolutionary Strategy (CMA-ES), an Approach to Investigate Differences in Channel Properties Between Neuron Subtypes. Frontiers in Neuroinformatics. 2018; 12. https://www.frontiersin.org/journals/neuroinformatics/articles/10.3389/fninf.2018.00047, doi: 10.3389/fninf.2018.00047.