# Bayesian-for-BOLD
This experiment will focus on developing simulations for blood oxygen level-dependent (BOLD) signals in functional Magnetic Resonance Imaging (fMRI) data that serve as the ground truth for Computational Parametric Mapping (CPM), a Bayesian approach for model-based fMRI experiments.

## Background 
<img align = "right" src="https://raw.githubusercontent.com/EstellaD/Bayesian-for-BOLD/master/img/introBOLD.png" width=700> [1]
<img align = "right" src="https://raw.githubusercontent.com/EstellaD/Bayesian-for-BOLD/master/img/ConvolveHRF.png" width=400> [2]

<img align = "right" src="https://raw.githubusercontent.com/EstellaD/Bayesian-for-BOLD/master/img/BOLDsignal.png" width=400> [3]


Give the shape of one single BOLD signal curve, we chose a conjugate prior for inverse gamma distribution. There are two steps for the Bayesian framework: 
1. First
2. Second

<img src="https://latex.codecogs.com/gif.latex?\alpha"/> is the power in steven's power law. We want to explore which value of alpha gives the best parameter recovery when designing experiments. 

<img src="https://latex.codecogs.com/gif.latex?\beta"/> This study is a simplified case with only one beta in the linear model. 



For more information, please refer to `img/SRP17_Dong_Conference_Post.pdf`.

## Running the simulation code

#### Requirements
This repository is implemented in MATLAB 3.7, with several dependency functions from [Statistical Parametric Mapping (SPM)] (https://www.fil.ion.ucl.ac.uk/spm/) package developed by UCL for neuroimging. 

#### Explanation of the file structure
`StarOneBeta.m` contains all visualization of the MAP of different power values. 


## References
[1] Arthurs, O. J., & Boniface, S. (2002). How well do we understand the neural origins of the fMRI BOLD signal? Trends in Neurosciences, 25(1), 27–31. https://doi.org/10.1016/S0166-2236(00)01995-0
[2] Zhang, L., Guindani, M., & Vannucci, M. (2015). Bayesian Models for fMRI Data Analysis. Wiley interdisciplinary reviews. Computational Statistics, 7(1), 21–41. https://doi.org/10.1002/wics.1339
[3] AD Elster, ELSTER LLC. (2020). BOLD and Brain Activity. Questions and Answers in MRI. http://mriquestions.com/does-boldbrain-activity.html
