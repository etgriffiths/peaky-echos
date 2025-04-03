# peaky-echos
Workflow for the detection and categorization of dolphin echolocation with peaks and notches present in their spectrum.

Emily T. Griffiths
emilytgrifffiths@ecos.au.dk
2025

The functions and example script presented can be used to replicate the decision tree and spectra template development used to discern and categorize white-beaked dolphins from other dolphins and noise in Skagerrak in 2021. Paper forthcoming. Parameters for template used in the paper included within the RDA file.

This repository also contains functions to further examine spectra, such as shoulder, which finds a prominent secondary peak or spectral lobe. This function is used in analysis outside of Skagerrak study (e.g. Griffiths et al., 2020), but in conjunction with other functions in this repository.

This workflow uses functions from the PamBinaries and PAMpal R packages to extract binary information per click, and ICES standards for measuring underwater noise.

https://github.com/TaikiSan21/PamBinaries

https://taikisan21.github.io/PAMpal/

https://github.com/ices-tools-prod/underwaternoise

Griffiths, E. T., Archer, F., Rankin, S., Keating, J. L., Keen, E., Barlow, J., & Moore, J. E. (2020). Detection and classification of narrow-band high frequency echolocation clicks from drifting recorders. The Journal of the Acoustical Society of America, 147(5), 3511-3522.
