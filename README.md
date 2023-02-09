# Electrochemical lithium-ion battery model and easy-to-apply parameterization procedures for fast-charge applications

This repository provides an electrochemical lithium-ion battery model for a commercially-available battery cell capable of simulating a battery system at various fast charging conditions
synchronously providing battery internal states (such as the cell's anode potential) for the investigation of fast charging control and optimization. A facile parameterization
procedure of the model has been derived, which can be adapted and transferred to other commercial lithium-ion battery cell's. All necessary tools and auxiliaries are provided.

The model is based on the Single Particle Model with electrolyte and thermal dynamics originally developed and published by the research group of Scott Moura, University of California, Berkeley, US,
updated and extended with a full parameter set, model validation for a commercial Sony US18650VTC5A cell (SiC/NCA) and parameter identification functions developed and published by

[Nikolaos Wassiliadis](mailto:nikolaos.wassiliadis@tum.de),<br/>
**[Institute of Automotive Technology](https://www.mos.ed.tum.de/mos/startseite/)**,<br/>
**[Technical University of Munich, Germany](https://www.tum.de/nc/en/)**,<br/>
https://doi.org/10.1016/j.est.2022.105951.

## Features
- In-depth analysis of half and full cell data of a commercially available lithium-ion battery (Sony US18650VTC5A).
- Novel hybrid parameter identification procedure with the required functions based on half cell and full cell measurements.
- Validation at charging currents up to 6C and various cell temperatures starting from -10°C up to 50°C.
- Validation from cell up to system level to demonstrate the model scalability.
- Sensitivity analysis of parameter variations on anode potential control.

## Usage of the battery model

We are very happy if you choose this model implementation and auxiliaries for your projects and provide all updates under GNU LESSER GENERAL PUBLIC LICENSE Version 3 (29 June 2007).
Please refer to the license file for any further questions about incorporating this battery model into your projects.
We are looking forward to hearing your feedback and kindly ask you to share bugfixes, improvements and updates on the parameterization or implementation.

## Requirements

The model was created with MATLAB 2019b and usage of the Curve Fitting and Optimization Toolbox. If you want to commit an updated version using another software release or a specific toolbox please give us a short heads-up. 

## How To Use

The main path contains the model as MATLAB implementation with the necessary parameter set. *spmet.m* is the main file to call the overall model.
All subfunctions relevant to the processing procedure of the main file can be found in *auxiliaries*.
The directory *parameter* contain the parameter set including the anode and cathode open-circuit potentials, the balancing and alignment and the kinetics fitting to the half and full cell.
The directory *measurement* serves as input folder for battery measurements.
All functions and measurements used for the parameterization procedure are provided in the directory *parameterization*.
The directory *sensitivity* provides all necessary files for a sensitivity analysis.

## Authors and Maintainers

- Nikolaos Wassiliadis, nikolaos.wassiliadis@tum.de
  - Idea and structure behind the model extensions and novel parameterization procedure.
  - Experimental investigation of the cell under study.
  - Supervision of the contributing student's theses.
  - Final revision and validation of the model.

## Contributions

The authors want to thank the following persons, who, with their brilliant contributions made this work possible.

- Leo Wildfeuer, leo.wildfeuer@tum.de
  - Two-electrode half cell measurements and analysis
- Manuel Ank, manuel.ank@tum.de
  - System validation measurement setup and analysis
- Andreas Bach, andreas.bach@tum.de
  - Implementation of the parameterization procedure as part of his master's thesis
- Matthias Wanzel, matthias.wanzel@tum.de
  - Parameter study and model implementation as part of his master's thesis
- Jakob Weiss, jakob.weiss@tum.de
  - Thermal parameterization during his master's studies
- Ann-Sophie Zollner, ann-sophie.zollner@tum.de
  - Thermal parameterization during her master's studies
  
## License

This project is licensed under the GPL License - see the LICENSE.md file for details.
The project is built upon the original contributions by the research group of Scott Moura, University of California, Berkeley, US,
which are themselves licensed under the GPL License. Files contain references and further information in the header.