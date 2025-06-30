# MELODY-Sand-Rubber-Mixture

This repository contains the simulation code and postprocessing tools used in the following studies on the micromechanical behavior of sand–rubber mixtures under shear:

- **Hu, Z., Shi, Y. H., Guo, N., & Yang, Z. X. (2022)**  
  *Micromechanical investigation of the shear behaviors of sand‒rubber mixtures using a multibody meshfree method.*  
  *Granular Matter, 24(3), 73.*  
  [DOI: 10.1007/s10035-022-01294-w](https://doi.org/10.1007/s10035-022-01294-w)

- **Guo, N., Shi, Y., Hu, Z., & Yang, Z. (2023)**  
  *Shear behaviour of sand–rubber mixtures: interpretation from energy transformation and force chain evolution.*  
  *Géotechnique, 75(1), 41–55.*  
  [DOI: 10.1680/jgeot.21.P.205](https://doi.org/10.1680/jgeot.21.P.205)

---

## 🔬 Overview

This code implements a **multibody meshfree method** for simulating the large-deformation behavior of sand–rubber mixtures under drained shear. The rubber particles are modeled as deformable hyperelastic bodies, while sand particles are represented as rigid bodies with frictional contact.

Key features:
- Rigid–deformable multibody coupling
- Implicit contact algorithm
- Energy tracking and force chain analysis
- Reproducible cases from published journal articles

---

## 📦 Dependencies
- MATLAB (for preprocessing scripts if needed)

---

## 🧰 Installation

For detailed instructions, please refer to [`MELODY_USER_GUIDE_3.94.pdf`](https://github.com/Yihao-Shi/MELODY-Sand-Rubber-Mixture/blob/main/MELODY_USER_GUIDE_3.94.pdf) in this repository.

## 📈 Postprocessing

Tools are provided in the `postprocessing/` folder for:

- Stress–strain curves
- Energy dissipation analysis
- Force chain visualization
- Particle motion and rotation tracking

---

## 🙏 Acknowledgments

This code is based on [MELODY 2D v3.94](http://guilhem.mollon.free.fr/Theme12_Eng.html) developed by **Guilhem Mollon** (INSA Lyon).  
We gratefully acknowledge the original author to provide the full user guide ([`MELODY_USER_GUIDE_3.94.pdf`](https://github.com/Yihao-Shi/MELODY-Sand-Rubber-Mixture/blob/main/MELODY_USER_GUIDE_3.94.pdf)) in this repository for reference.

---

## 📄 License

This project is released under the **MIT License** unless stated otherwise in individual source files.

---

## 📬 Contact

**Yihao Shi**  
📧 shiyh99 [at] gmail.com 
🔗 [https://github.com/Yihao-Shi](https://github.com/Yihao-Shi)
