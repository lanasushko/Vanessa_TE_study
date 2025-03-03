# Code repository for annotation and analysis of TEs in Vanessa butterflies

This repository contains code for TE annotation of 4 buttefly species and the subsequent analyses that have been performed.

### 01_TE_discovery_and_library_construction

De novo TE libraries for the four species were constructed independently using RepeatModeler2, version 2.0.4 for V. cardui and V. atalanta; and version 2.0.1 for V. tameamea and A. io with LTR structural module as in `repeatmodeler.sh`. MITE families were detected and added separately using MITE Tracker pipeline [`mitetracker.sh`]. MCHelper was run in automatic mode on the four TE libraries to extend fragmented consensus sequences and remove false positives based on structural criteria [`mchelper.sh`].


