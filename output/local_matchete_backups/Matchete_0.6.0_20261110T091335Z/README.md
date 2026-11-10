# Matchete

[![arXiv](https://img.shields.io/badge/arXiv-2212.04510-00aa00.svg)](https://arxiv.org/abs/2212.04510)

Matchete is a [Mathematica](https://www.wolfram.com/mathematica/resources/) package aimed at facilitating the functional matching procedure for generic weakly coupled UV models with a mass power counting. It is the first package fully automating functional one-loop matching computations. It is built upon and superseeds the [SuperTracer](https://gitlab.com/supertracer/supertracer) package [\[arXiv:2012.08506\]](https://arxiv.org/abs/2012.08506), which is used for the evaluation of the functional supertraces.

Matchete provides a simple and user-friendly interface for entering the Lagrangian of a generic UV theory. First of all, the user has to specify all the (gauge) symmetry groups and representations of the model, then all fields and couplings can be defined. Afterwards, the Lagrangian can be written very close to a pen-and-paper form.

The `Match` routine can then be applied to the UV Lagrangian to obtain the corresponding Effective Field Theory, where all heavy degrees of freedom have been integrated out, either at tree level, or at one loop.

The resulting EFT Lagrangian contains a (very) large number of redundant operators. These can be removed automatically by applying the routines `GreensSimplify`, which reduces the output to an off-shell Green's basis, and `EOMSimplify`, which performs field redefinitions. The fully simplified result will then be in a (near-)basis of the corresponding EFT operator space.

More details can be found at [matchete.gitlab.io](https://matchete.gitlab.io/) and a [model repository](https://gitlab.com/matchete/model-database) can be consulted for inspiration and convenience.

---

## Reference

If you use Matchete please cite: [\[arXiv:2212.04510\]](https://arxiv.org/abs/2212.04510).

---

## Installing and loading of the package

The simplest way to download and install Matchete is to run the following command in a Mathematica session:

> Import["https://gitlab.com/matchete/matchete/-/raw/master/install.m"]

This will download and install Matchete in the Applications folder of Mathematica's base directory. To load Matchete use the command:

> <<Matchete\`

The complete set of routines and usage examples can be found in the ancillary documentation notebooks. This repository also contains several example notebooks demonstrating the usage of Matchete for the one-loop matching of the following theories:
* Vector-like fermion toy model
* Singlet scalar extension of the SM
* Vector-like charged lepton singlet extension of the SM

---

## Authors

* **Javier Fuentes-Martín** - *Universidad de Granada*
* **Matthias König** - *University of Mainz*
* **Julie Pagès** - *SLAC National Accelerator Laboratory*
* **Anders Eller Thomsen** - *University of Bern*
* **Felix Wilsch** - *RWTH Aachen University*

## Contributors
* **Lukas Born** - *University of Bern*
* **Luis Hourtz** - *Paul Scherrer Institute (PSI)* 
* **Tamilarasan Ketheeswaran** - *RWTH Aachen University*
* **Javier Lopéz-Miras** - *Universidad de Granada*
* **Adrián Moreno-Sánchez** - *Universidad de Granada*
* **Philip Weber** - *RWTH Aachen University*

---

## Bugs and feature requests

Please submit bugs and feature requests using GitLab's [issue system](https://gitlab.com/matchete/matchete/-/issues).

---

## License

MATCHETE is free software under the terms of the GNU General Public License v3.0.

---

## Acknowledgments

We thank José Santiago for his help with the cross-checks of the vector-like lepton example using matchmakereft, and Mikael Chala and Guilherme Guedes for their help in validating the `Matchotter` module.
