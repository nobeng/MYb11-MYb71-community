This repository contains the code and data for the project and publication "Role of microbial life history strategy in shaping the characteristics and evolution of host-microbiota interactions" (Obeng et al. 2025). A corresponding release has been published at Zenodo. All analyses performed in R, their corresponding data and output from statistical analyses are archived in this repository.

## Table of Contents
•	Abstract
•	Repository Structure
•	Getting Started
•	Reproducing Results
•	License
•	Contact
•	Citation

## Abstract

Many host-associated microbes are transmitted between individual hosts via the environment and, therefore, need to succeed both within a host and a connected environmental habitat. These microbes might invest differentially into the two habitats, potentially leading to fitness trade-offs and distinct life history strategies that ultimately shape the host-associated microbial communities. In this study, we investigated how the presence of distinct bacterial life history strategies affects microbiota characteristics along a host-associated life cycle, using the nematode host *Caenorhabditis elegans* and two naturally associated bacteria, *Pseudomonas lurida* and *Ochrobactrum vermis*, as an experimentally tractable model. Based on genomic life history prediction and experimental fitness characterizations, we identified distinct ecological strategies for the bacteria: while *P. lurida* dominated the free-living environment, *O. vermis* was more abundant in the host. Using mathematical modelling, experimental evolution, and whole genome sequencing, we next assessed whether the two distinct ecological strategies influence further adaptation to the host-associated life cycle. We found that (i) the host-specialist *O. vermis* did not further adapt to the two habitats, whereas (ii) the initially better environmental competitor *P. lurida* did adapt to the life cycle, leading to its increased abundance in both environment and host. Evolutionary adaptation of *P. lurida* caused a shift in microbiota composition in the host, which in turn, resulted in a significant increase in host fitness. Overall, our results highlight the role of microbial life history strategies in shaping the characteristics and evolution of host-microbe interactions and point to a potential selective advantage of better environmental competitors.

Authors: Nancy Obeng, Johannes Zimmermann, Anna Czerwinski, Janina Fuß, Hinrich Schulenburg
Acknowledgements: We thank Florence Bansept, Brendan Bohannan, Peter Deines, and the Schulenburg lab for critical feedback, Anna Czerwinski for lab support, the Kiel BiMo/LMB for access to their core facilities. Funding was provided by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation), Project-ID 261376515 – SFB 1182, Projects A4.3 and B4.3 (NO, HS), the DFG Research Infrastructure NGS_CC project 407495230 (JF) as part of the Next Generation Sequencing Competence Network project 423957469, the International Max-Planck Research School for Evolutionary Biology (NO), and the Max-Planck Society (Fellowship to HS).

The repository includes scripts and data necessary to reproduce the analyses and figures presented in the associated publication.

## Repository Structure
R scripts are saved under the main branch. Most input data is saved under the main branch, alternatively in "data". Output from statistical analyses can be found under "stats", output for life history predictions under "plots".

## Getting Started
1.	Clone the repository:
bash
git clone https://github.com/nobeng/MYb11-MYb71-community.git
cd MYb11-MYb71-community

2.	Install R and required packages:
•	R version used before upload: v4.4.1
•	All required R packages are listed at the top of each script. Install them using: install.packages(c("package1", "package2", ...))

3.	Run the analysis
•	Execute the scripts

## Reproducing Results
•	All scripts are self-contained and require no special instructions beyond installing dependencies.
•	Input data files are provided in the repository.

## License
This project is licensed under the Creative Commons Attribution 4.0 International (CC BY 4.0) license.

## Contact
For questions or further information, please raise a GitHub Issue or Discussion.

## Citation
If you use this code or data, please cite the associated publication: [Full citation to be added upon publication].
