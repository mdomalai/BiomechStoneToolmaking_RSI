This repository contains the scripts associated with the study "Biomechanical demands of percussive techniques in the context of early stone toolmaking" (Macchi R., Daver G., Brenet  M., Prat S., Hugheville L., Harmand S., Lewis J., Domalain M.) published in the Journal of the Royal Society Interface, 2021. https://doi.org/10.1098/rsif.2020.1044

The scripts were originally developped to process the data used in this study and available on datadryad.org (https://doi.org/10.5061/dryad.cfxpnvx51).

To start: Put the EMG, Kinematic, RMC1_marker and EMG_maximal_voluntary_contraction files in different folders (there should result in 4 folders plus 1 for the backup files).

You can now run the first code that processes all the EMG data (detects the cycle phases, performs EMG envelopes in addition to iEMG values for each muscle).
During the processing, you have to choose the pathway where the data are located:
First -> the pathway of the EMG data
Second -> the pathway of R_MC1 markers 
Third -> the patheway of maximal vouluntary contraction EMG 
Fourth -> the patheway of the backup files

Then, you can run the second code that processes the kinematic data (i.e computes the normalized cycles). 

Each code plots a visualization of each muscle and kinematic cyle that allows to manually remove cycle errors.

May you have any question, don't hesitate to contact : robin.macchi@gmail.com


![Image of Macchi et al. RSI 2021](https://github.com/mdomalai/BiomechStoneToolmaking_RSI/blob/main/ImageGitHub.png)











