Stochastic model of ryanodine type 2 and 3 receptor activation by calcium and gating by calmodulin
=======================


Introduction
-----------------------
This repo contains two kinetic models: one of the ryanodine type 2 receptor (RyR2) inactivation by calmodulin and activation by calcium and ryanodine type 3 receptor (RyR3) activation by calcium and calmodulin. RyR2 is expressed mostly in cardiac myocites and has been quite thoroughly investigated. However, there is much less known about the RyR3, which is express throught the body. We used the known measurements of affinity of calmodulin binding to CaMBD2 domain of RyR2 to built both RyR2-CaM and RyR3-CaM [1, 2, 3] and dissociation rates [4, 5], because CaMBD2 is highly conserved for both ryanodine receptor types [6].

Models have been fitted to experimental data. Calcium activation of RyR2 unbound to calmodulin was fitted to RyR2 opening time course from [7], and verfied by evaluating the fit to open probability, mean open and closed time provided by [3]. Calcium activation of CaM-bound RyR2 open probability was fitted to appropriate data from [3] and the mean open and closed time was validated to best match data from [3]. 

Models are implemented using NeuroRD, a stochastic reaction-diffusion simulator implementing Gillespie algorith with asynchronous tau-leaping.


[1] Lau K, Chan MMY, Van Petegem F. Lobe-Speciﬁc Calmodulin Binding to Different Ryanodine Receptor Isoforms. Biochemistry. 2014; 53(5):932–946. https://doi.org/10.1021/bi401502x, doi: 10.1021/bi401502x, pMID:
24447242.

[2] Søndergaard MT, Liu Y, Guo W, Wei J, Wang R, Brohus M, Overgaard MT, Chen SRW. Role of cardiac ryanodine receptor calmodulin-binding domains in mediating the action of arrhythmogenic calmodulin N-domain mutation N54I. The FEBS Journal. 2020; 287(11):2256–2280. https://febs.onlinelibrary.wiley.com/doi/abs/10.1111/febs.15147, doi: https://doi.org/10.1111/febs.15147.

[3] Xu L, Meissner G. Mechanism of Calmodulin Inhibition of Cardiac Sarcoplasmic Reticulum Ca2+ Release Channel (Ryanodine Receptor). Biophysical Journal. 2004; 86(2):797–804. https://www.sciencedirect.com/science/article/pii/S0006349504741557, doi: https://doi.org/10.1016/S0006-3495(04)74155-7.

[4] Balshaw DM, Xu L, Yamaguchi N, Pasek DA, Meissner G. Calmodulin Binding and Inhibition of Cardiac Muscle Calcium Release Channel (Ryanodine Receptor)*.
Journal of Biological Chemistry. 2001; 276(23):20144–20153. https://www.sciencedirect.com/science/article/pii/S0021925819404547, doi:https://doi.org/10.1074/jbc.M010771200.

[5] Wu X, Bers DM. Free and bound intracellular calmodulin measurements in cardiac myocytes. Cell Calcium. 2007; 41(4):353–364. https://www.sciencedirect.com/science/article/pii/S0143416006001618, doi:https://doi.org/10.1016/j.ceca.2006.07.011.

[6] Brohus M, Søndergaard MT, Wayne Chen SR, van Petegem F, Overgaard MT. Ca2+-dependent calmodulin binding to cardiac ryanodine receptor (RyR2) calmodulin-binding domains. Biochemical Journal. 2019 01; 476(2):193–209. https://doi.org/10.1042/BCJ20180545, doi: 10.1042/BCJ20180545.

[7] Zahradníková A, Zahradník I, Györke I, Györke S. Rapid Activation of the Cardiac Ryanodine Receptor by Sub-
millisecond Calcium Stimuli. Journal of General Physiology. 1999 11; 114(6):787–798. https://doi.org/10.1085/
jgp.114.6.787, doi: 10.1085/jgp.114.6.787.