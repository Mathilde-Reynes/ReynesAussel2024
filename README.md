# [Re] Model of Thalamocortical Slow-Wave Sleep Oscillations and Transitions to Activated States

**Authors:** Mathilde Reynes, Amélie Aussel  
**Contact:** mathilde.reynes@u-bordeaux.fr

## 1. References
The model used in this simulator is detailed in an article [1] that is yet to be published. The simulator relies on the Brian2 libraries for Python [2].


## 2. Requirements
This simulator was developed using Python 3.8 and the following packages: brian2=2.5.1, python=3.8.19, numpy=1.24.4. It was developed on Windows 11 and tested on Linux. The specifications for the Windows environment are available in the environment.yml file in the Model folder of this GitHub repository. You can create a similar virtual environment with Conda using the command:

```bash
conda env create -f environment.yml
```

This simulator was tested on the latest python (3.11.7) and brian2 (2.7) versions and reproducibility was ensured.

## 3. User Interface
The model is organized with each component in its own Python file:

### Cells:
- Soma_eqs.py: Equations for pyramidal and inhibitory interneuron axosomatic compartments.
- Dendritic_eqs.py: Equations for pyramidal and inhibitory interneuron dendritic compartments.
- RE_eqs.py: Equations for reticular cells.
- TC_eqs.py: Equations for thalamic relay cells.

### Synapses:
- Synapses.py: Equations for AMPA, NMDA, GABA_A, and GABA_B synapses.

### Thalamocortical sub-systems:
- Cortical_layer.py: Defines the cortical compartment, including pyramidal and interneurons, and their synapses.
- Thalamus.py: Defines the thalamic compartment, including thalamic relay and reticular cells, and their synapses.

### Simulation:
- Thalamo_cortical.py: Simulates the full model through the function thalamocortical_network(seed_value,analyze_speed,fig_number,raw_data,plot_figure) with: 
    - analyze_speed (True/False) allows the user to compute the mean and std propagation speed of cortical up-states
    - fig_number allows the user to run the simulation with the parameters necessary to reproduce a specific figure from the article; possible values are "5, 7, 8, 9-A1, 9-A2, 9-A3, 9-B1, 9-B2, 9-B3, 10-1, 10-2, 10-3, 11-1, 11-2, 12-A1, 12-A2, 12-B1, 12-B2, 13-A, 13-B, 14, 15-A, 15-B, 15-C, 15-D, 16-A1, 16-A2, 16-A3, 16-A4, 16-B1, 16-B2, 16-B3, 16-B5, 16-B4, 16-B6, 16-B7, 17-A1, 17-A2, 17-A3, 18-A1, 18-A2, 18-A3, 17-B1, 17-B2, 17-B3, 18-B1, 18-B2, 18-B3, 19, S1, S2, S3-A1, S3-A2, S3-A3, S3-A4, S3-B1, S3-B2, S3-B3, S3-B4, S3-B5, S3-B6, S3-B7"
    - raw_data (True/False) allows the user to select if raw data (.txt files) should be saved
    - plot_figure (True/False) allows the user to plot the figure corresponding to the fig_number

### Files with functions that manage parameter updates, save raw data, and generate plots for the figures : 
- Figure_conditions.py : Designed to configure the model with specific parameters adjusted for different figures. 
- Figure_plotting.py : Designed to plot a specific figure of the article.
- Figure_rawdata.py : Designed to save the raw data necessary for a specific figure.

## 4. Additional Files
Code for plotting figures 2, 3, 4 using Bazhenov et al. (2002) original data is available in the ‘Figures’ folder on the GitHub repository. [Dataset 1](https://doi.org/10.5281/zenodo.13308394) generated from the present model and [Dataset 2](https://doi.org/10.5281/zenodo.13308394) generated from the original model are both accessible via Zenodo. After download, both folders should be placed in \Reynes_Aussel, alongside Figures and Model folders.


References

[1] Reynes, M., & Aussel, A., to be published

[2] Stimberg, M., Brette, R., & Goodman, D. F. M. (2019). “Brian 2, an Intuitive and Efficient Neural Simulator.” eLife, 8, e47314. doi: 10.7554/eLife.47314