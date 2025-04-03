<!-- PROJECT LOGO -->
<br />
<p align="center">
  <a href= "">
  <h3 align="center"> Contact percolation governs collecive motility via mechanochemical feedback in eterogeneouse breast cancer
</h3>
  </a>
  <p align="center">
  Leonardo Barzaghi<sup>* 1,2 </sup> , Chiara Guidolin<sup>3</sup>, Camillo Mazzella<sup>1</sup>, Joseph Ackermann<sup>4</sup>, Edoardo Bellini<sup>1,2</sup>, Stefano Villa<sup>3,5</sup>, Emanuela Frittoli<sup>1</sup>, Brenda Green<sup>1</sup>, Andrea Palamidessi<sup>1</sup>, Angela Cattaneo<sup>6</sup>, Angela Bachi<sup>1</sup>, Valeria Cancila<sup>7</sup>, Claudio Tripodo<sup>1,2,7</sup>, Simona Polo<sup>1,2</sup>, Roberto Cerbino<sup>8</sup>, Raphaël Voituriez<sup>4</sup>, Fabio Giavazzi </sup>*3, and Giorgio Scita <sup>†1, 2</sup>
  <details open="open">
  <summary>Reference</summary>
  <ol>
  <li>
      <a href="1"> IFOM ETS, the AIRC Institute of Molecular Oncology, Milan, Italy </a>
  </li>
    <li>
      <a href="2"> Department of Oncology and Haemato-Oncology, University of Milan, Milan, Italy </a>
    </li>
    <li>
      <a href="3"> Department of Medical Biotechnology and Translational Medicine, University of Milan, Segrate, Italy </a>
    </li>
    <li>
      <a href="4"> Laboratoire Jean Perrin, CNRS/Sorbonne Université, Paris, France </a>
    </li>
    <li>
      <a href="5"> Max Plank Institute for Dynamics and Self-Organization, Göttingen, Germany </a>
    </li>
    <li>
      <a href="6"> Cogentech SRL Benefit Corporation, Milan, Italy </a>
    </li>
    <li>
      <a href="7"> Department of Health Sciences, Human Pathology Section, University of Palermo School of Medicine, Palermo, Italy </a>
    </li>
    <li>
      <a href="8"> Faculty of Physics, University of Vienna, Vienna, Austria </a>
    </li>
  </ol>
</details>
  </p>
</p>
<p align="center">
</p>
<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary>Table of Contents</summary>
  <ol>
  <li>
      <a href="#About the paper">About the paper</a>
  </li>
    <li>
      <a href="#Pipelines">Pipelines</a>
    </li>
    <li>
      <a href="#Contacts">Contacts</a>
    </li>
  </ol>
</details>

<!-- About the paper -->
## About the paper

Breast cancer progression is driven by dynamic changes in tissue mechanics that promote phase transitions toward collective cell movement, immune activation, and invasion. While previous investigations primarily employed model systems of genetically identical cells, breast carcinoma consists of heterogeneous cell populations varying in genetic and mechanical traits. The complex interplay between this diversity and collective cell dynamics, affecting phase transitions and tumor progression, remains largely unexplored.Here, we interpreted our experimental data through contact percolation theory to uncover a chemo-mechanical switch underlying a phase transition via the emergence of flocking motility in heterogeneous breast cancer tissues. Using engineered mixtures of motile RAB5A-expressing cells and immotile controls, we demonstrate that interconnected system-spanning clusters of fluidized RAB5A cells induce a phenotypic switch in neighboring cells, enabling them to acquire full-scale, flocking-like motility and reprogram their transcriptional
state. This transition is accompanied by the activation of pro-inflammatory gene programs. Our theoretical modeling supports a mechanism by which local cell composition drives a motility switch critical for the emergence of coordinated, system-wide movement. These findings highlight the role of mechanical heterogeneity in tumor progression and identify contact percolation as a fundamental mechanism driving dynamic tissue reorganization and immune modulation at the root of breast cancer invasion and metastasis.

<!-- Pipelines -->
## Pipelines
### RNA_seq analisys
In this [folder](./Pipelines/RNA_seq/) are located the pipelines, bash and R script used to analyze all the RNA-seq data [deposited]() . Within the folder are present:

* `.ipynb` file, custom R script used to analyze a specific dataset.

* `.sh` files, bash script used to launch specific softwares and [nf-core/RNA-seq](https://github.com/nf-core/rnaseq/tree/3.14.0) pipeline.

* `.csv` files, text file containing metadata information of the samples.

### Image analisys
In this [folder](./Pipelines/Image_analisys/) are present the pipelines, python scripts used to analyze immunofluorescence images for cell segmentation. Within the folder are present:

* `launcher_segmentation.py` python script that launches the segmentation analisys according to user-defined parameters in the file [parametri.yaml](./Pipelines/Image_analisys/parametri.yaml) 

* `segmentation.py` python script that segments cells and nuclei using [cellpose](https://www.nature.com/articles/s41592-020-01018-x.epdf?sharing_token=yrCA1mB-y9TR8-RC8w_CPdRgN0jAjWel9jnR3ZoTv0Ms-A3TbUG5N7s_6d3I7lMImMDE6cyl-17ubiknffX50r-dX1un0XSIQ2PGYWsCV1du16fIaipcHNxste8FMByEL75Ek_S2_UEVkSk7lCFllWEVogGWJwmQkBC9uKq9UEA%3D) and extracts general shape parameters.

* `mixedpop_segmentation.py` python script that segments cells and nuclei from mixed population of cells with different nuclear florescent markers and extracts general shape parameters.

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
