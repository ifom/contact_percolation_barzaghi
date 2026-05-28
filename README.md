<!-- PROJECT LOGO -->
<br />
<p align="center">
  <a href="https://www.researchsquare.com/article/rs-6057847/v1" target="_blank">
    <strong style="font-size: 1.2em;">
      Contact percolation governs collective motility via mechanochemical feedback in heterogeneous breast cancer
    </strong>
  </a>
</p>
  <p align="center">
    Leonardo Barzaghi<sup>†,1,2</sup>, Camillo Mazzella<sup>†,1</sup>, Chiara Guidolin<sup>3</sup>,  Joseph Ackermann<sup>4</sup>, Edoardo Bellini<sup>1,2</sup>,  Andrew E. Massey<sup>5</sup>, Stefano Villa<sup>3,6</sup>, Emanuela Frittoli<sup>1</sup>, Brenda Green<sup>1</sup>, Andrea Palamidessi<sup>1</sup>, Angela Cattaneo<sup>7</sup>, Angela Bachi<sup>1</sup>, Valeria Cancila<sup>8</sup>, Claudio Tripodo<sup>1,2,8</sup>, Simona Polo<sup>1,2</sup>, Roberto Cerbino<sup>9</sup>, Alexander X. Cartagena-Rivera<sup>5</sup>, Raphaël Voituriez<sup>4</sup>, Fabio Giavazzi<sup>#3</sup>, and Giorgio Scita<sup>#1,2</sup>
  </p>
</p>

† These authors contributed equally.
    # Correspondence should be addressed to Fabio Giavazzi (fabio.giavazzi@unimi.it) or Giorgio Scita (giorgio.scita@ifom.eu).


<details open="open">
  <summary>Reference</summary>
  <ol>
    <li><a href="#">IFOM ETS, the AIRC Institute of Molecular Oncology, Milan, Italy</a></li>
    <li><a href="#">Department of Oncology and Haemato-Oncology, University of Milan, Milan, Italy</a></li>
    <li><a href="#">Department of Medical Biotechnology and Translational Medicine, University of Milan, Segrate, Italy</a></li>
    <li><a href="#">Laboratoire Jean Perrin, CNRS/Sorbonne Université, Paris, France</a></li>
    <li><a href="#">Section on Mechanobiology, National Institute of Biomedical Imaging and Bioengineering, National Institutes of Health, Bethesda, MD, USA</a></li>
    <li><a href="#">Max Planck Institute for Dynamics and Self-Organization, Göttingen, Germany</a></li>
    <li><a href="#">Cogentech SRL Benefit Corporation, Milan, Italy</a></li>
    <li><a href="#">Department of Health Sciences, Human Pathology Section, University of Palermo School of Medicine, Palermo, Italy</a></li>
    <li><a href="#">Faculty of Physics, University of Vienna, Vienna, Austria</a></li>

  </ol>
</details>

<details open="open">
  <summary>Table of Contents</summary>
  <ol>
    <li><a href="#about-the-paper">About the paper</a></li>
    <li><a href="#pipelines">Pipelines</a></li>
    <li><a href="#contacts">Contacts</a></li>
  </ol>
</details>

<!-- About the paper -->
## About the paper

Breast cancer progression is driven by dynamic changes in tissue mechanics that promote phase transitions toward collective cell movement, immune activation, and invasion. Whilethese processes have been mainly investigated in model systems of genetically identical cells, breast carcinomas are intrinsically heterogeneous, comprising cell populations with distinct genetic and mechanical properties. How this heterogeneity shapes collective dynamics and tumor progression remains poorly understood. Here, we interpret our experimental data through contact percolation theory to uncover a mechano-chemical switch underlying the emergence of collective flocking motility in heterogeneous breast cancer tissues. Using engineered mixtures of motile RAB5A-expressing cells and immotile control cells, we show that system-spanning clusters of fluidized RAB5A cells induce an abrupt morphodynamic switch in neighboring control cells, enabling them to acquire polarized protrusions and full-scale flocking-like motility. Theoretical modeling supports a mechanism by which local cell composition governs this transition and the onset of coordinated, system-wide motion, in a process that depends on RAC1-driven protrusive activity. Mechanistically, control cells in mixed monolayers above the percolation threshold exhibit activation of an E- cadherin-dependent EGFR–MAPK signaling axis, whereas E-cadherin silencing suppresses MAPK activation, impairs protrusion alignment, and disrupts flocking. This transition is accompanied by collective mechanical remodeling, as control cells progressively soften and become more elongated as the fraction of neighboring RAB5A cells increases, consistent with enhanced tissue fluidization. In parallel, control cells activate a STAT1/STAT2-dependent inflammatory program, fostered by the combined action of paracrine signaling and contact-dependent amplification. Together, these findings support a model in which contact percolation and mechano-chemical coupling contribute to dynamic tissue reorganization, inflammatory reprogramming, and invasive behavior in heterogeneous breast cancer.

<!-- Pipelines -->
## Pipelines
### RNA_seq analisys
In this [folder](./Pipelines/RNA_seq/) are located the pipelines, bash and R script used to analyze all the RNA-seq data [deposited]() . Within the folder are present:

* `.ipynb` file, custom R script used to analyze a specific dataset.

* `.sh` files, bash script used to launch specific softwares and [nf-core/RNA-seq](https://github.com/nf-core/rnaseq/tree/3.14.0) pipeline.

* `.csv` files, text file containing metadata information of the samples.

* `.yaml` file, cotaining all requirments to build a conda enviroment to reproduce the analisys and plots of the RNA_seq.

### CellposeSAM_Segmentation_Pipeline

This [folder](./Pipelines/CellposeSAM_Segmentation_Pipeline/) contains Python scripts and pipelines for analyzing immunofluorescence images and performing cell segmentation. For detailed information regarding the pipeline functionality, see the [README.md](./Pipelines/CellposeSAM_Segmentation_Pipeline/README.md).

### Analisys of cell Motility
In this [folder](./Pipelines/analisys_of_cell_motility) are present the pipelines, and MATLAB scripts used to analyze maps of the flow field in the cell monolayer are obtained by analyzing phase-contrast movies with PIV. Within the folder are present :

* [`README.md`](./Pipelines/analisys_of_cell_motility/README.md) file contain all instruction to run the pipeline regarding cell motility

The overall image analysis is divided in three main blocks:

1. **PIV**: preliminary analysis to identify the peak of maximum motility of 100% RAB5A population

2. **STATIX**: nuclei segmentation and characterization of the cell pattern from the segmented cell nuclei

3. **DYNAMIX**: single-cell tracking and quantification of the main dynamical parameters of interest from cell trajectories



<!-- CONTACTS -->
## Contacts

Leonardo Barzaghi - <leonardo.barzaghi@ifom.eu>

Edoardo Bellini - <edoardo.bellini@ifom.eu>

Giorgio Scita - <giorgio.scita@ifom.eu>

<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/othneildrew/Best-README-Template.svg?style=for-the-badge
[contributors-url]: https://github.com/othneildrew/Best-README-Template/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/othneildrew/Best-README-Template.svg?style=for-the-badge
[forks-url]: https://github.com/othneildrew/Best-README-Template/network/members
[stars-shield]: https://img.shields.io/github/stars/othneildrew/Best-README-Template.svg?style=for-the-badge
[stars-url]: https://github.com/othneildrew/Best-README-Template/stargazers
[issues-shield]: https://img.shields.io/github/issues/othneildrew/Best-README-Template.svg?style=for-the-badge
[issues-url]: https://github.com/othneildrew/Best-README-Template/issues
[license-shield]: https://img.shields.io/github/license/othneildrew/Best-README-Template.svg?style=for-the-badge
[license-url]: https://github.com/othneildrew/Best-README-Template/blob/master/LICENSE.txt
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://linkedin.com/in/othneildrew
