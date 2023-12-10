Code repository for our version of ["Total Variation Mode Decomposition"](https://link.springer.com/chapter/10.1007/978-3-030-75549-2_5) by [Ido Cohen](https://idocohen.co.il/), [Tom Berkov](https://www.vision-and-sensing.com/copy-of-jonathan-brokman), [Guy Gilboa](https://guygilboa.net.technion.ac.il/), at Technion - Israel Institute of Technology. [Download](https://github.com/IdoCohen5743/Total-Variation-Mode-Decomposition/archive/refs/heads/main.zip)

In this work we analyze the Total Variation (TV) flow applied to one dimensional signals. We formulate a relation between Dynamic
Mode Decomposition (DMD), a dimensionality reduction method based on the Koopman operator, and the spectral TV decomposition. DMD is
adapted by time rescaling to fit linearly decaying processes, such as the TV flow. For the flow with finite subgradient transitions, a closed form 10 solution of the rescaled DMD is formulated. In addition, a solution to the TV-flow is presented, which relies only on the initial condition and its corresponding subgradient. A very fast numerical algorithm is obtained which solves the entire flow by elementary subgradient updates.

<p align="center"><img src="https://i.imgur.com/7VbdaoN.png" | height=250></p>
TV flow decomposed into nonlinear modes (using time-rescaled-DMD):

<p align="center"><img src="https://i.imgur.com/UKDQt6Z.png" | height=140 ></p>

Comparison of original TV flow vs. Time-Rescaled TV flow:
<p align="center"><img src="https://i.imgur.com/6Q5hDU4.gif"></p>




## Requirements
This code has been implemented using Matlab-2018b

## Data & directory tree

```
code_submission 	# code and data together
├── fast_experiment1D_final.m  		 # Full experiment script for toy example (Three pulses)
├── fast_experiment2D_zebra_final.m  # Full experiment script for harder example (Line from zebra image)
├── proj_tvl2.m          			 # Helper function for ss_freq_tv_evolve
├── ss_freq_tv_evolve.m   			 # Mathematical Engine of TV - evolve TV flow
├── TV_data1D.mat         		     # Data for toy example (three pulse) experiment
├── TV_ZebraLine.mat      			 # Data for zebra experiment
├── zebra_media_gmu.jpg   			 # Image from which the input signal for zebra experiment was taken
├── Additional Data       			 # More Data can be added here
```

## Run algorithm
```bash
# Run toy example experiment
fast_experiment1D_final
```

```bash

# For fast gradient calculation run
subgradient.m
```

```bash
# Run zebra row experiment
fast_experiment2D_zebra_final
```
Full computational mode:
```bash
# In order to run either code without relying on pre-saved data in .mat files, please change the loadData flag in line 3 to 0
loadData = 0;
```

## Runtime
Running in full computational mode (without loading data) might change from machine to machine. Results presented for running experiment on a 8-Gen. Core i7 laptop, with 16BG of RAM.

- A full run of the basic toy example experiment takes ~8.5 Minutes
- A full run of the advanced experiment with the row from the zebra image takes ~2 Hours

The above mentioned long times are due to usage of standard method for computing subgradients iteratively, which is the method we compare to.

- All other run modes take less than 1 minute



## Citation
If you find our work useful, please cite our paper:
```bash
@InProceedings{10.1007/978-3-030-75549-2_5,
author="Cohen, Ido
and Berkov, Tom
and Gilboa, Guy",
editor="Elmoataz, Abderrahim
and Fadili, Jalal
and Qu{\'e}au, Yvain
and Rabin, Julien
and Simon, Lo{\"i}c",
title="Total-Variation Mode Decomposition",
booktitle="Scale Space and Variational Methods in Computer Vision",
year="2021",
publisher="Springer International Publishing",
address="Cham",
pages="52--64",
isbn="978-3-030-75549-2"
}
and

@article{cohen2021total,
  title={Total-Variation Mode Decomposition},
  author={Cohen, Ido and Berkov, Tom and Gilboa, Guy},
  journal={arXiv preprint arXiv:2105.10044},
  year={2021}
}

```

Feel free to place issues here or contact me via the e-mail in my personal page.
