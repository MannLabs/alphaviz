![Pip installation](https://github.com/MannLabs/alphaviz/workflows/Default%20installation%20and%20tests/badge.svg)
![GUI and PyPi releases](https://github.com/MannLabs/alphaviz/workflows/Publish%20on%20PyPi%20and%20release%20on%20GitHub/badge.svg)
[![Downloads](https://pepy.tech/badge/alphaviz)](https://pepy.tech/project/alphaviz)
[![Downloads](https://pepy.tech/badge/alphaviz/month)](https://pepy.tech/project/alphaviz)
[![Downloads](https://pepy.tech/badge/alphaviz/week)](https://pepy.tech/project/alphaviz)
[![Documentation Status](https://readthedocs.org/projects/alphaviz/badge/?version=latest)](https://alphaviz.readthedocs.io/en/latest/?badge=latest)


# AlphaViz
**AlphaViz** is a cutting-edge browser-based interactive visualization tool allowing to visualize the processed mass spectrometry data acquired with **Bruker** instrument. The **AlphaViz** dashboard facilitates easy quality control of your analyzed samples and a clear inspection of the raw data of significant peptides/proteins.

To enable all hyperlinks in this document, please view it at [GitHub](https://github.com/MannLabs/alphaviz).

* [**About**](#about)
* [**License**](#license)
* [**Installation**](#installation)
  * [**One-click GUI**](#one-click-gui)
  * [**Pip installer**](#pip)
  * [**Developer installer**](#developer)
* [**Usage**](#usage)
  * [**GUI**](#gui)
  * [**CLI**](#cli)
  * [**Python and jupyter notebooks**](#python-and-jupyter-notebooks)
* [**Troubleshooting**](#troubleshooting)
* [**Citations**](#citations)
* [**How to contribute**](#how-to-contribute)
* [**Changelog**](#changelog)

---
## About

Software tools such as **MaxQuant** or **DIA-NN** identify and quantify high amounts of proteins. After downstream processing in **Perseus**, **MSstats** or the **Clinical Knowledge Graph**, differentially expressed proteins become possible candidates for biomarker discovery. **AlphaViz** is an automated visualization pipeline to link these identifications with the original raw data and easily assess their individual quality or the overall quality whole samples.

An open-source Python package of the AlphaPept ecosystem from the [Mann Labs at the Max Planck Institute of Biochemistry](https://www.biochem.mpg.de/mann). This project is built purely in Python using a new cutting-edge [Holoviz ecosystem](https://holoviz.org/index.html) and Plotly library to create interactive dashboards and plots.

---
## License

AlphaViz was developed by the [Mann Labs at the Max Planck Institute of Biochemistry](https://www.biochem.mpg.de/mann) and is freely available with an [Apache License](LICENSE.txt). External Python packages (available in the [requirements](requirements) folder) have their own licenses, which can be consulted on their respective websites.

---
## Installation

AlphaViz can be installed and used on all major operating systems (Windows, macOS, Linux).
There are three different types of installation possible:

* [**One-click GUI installer:**](#one-click-gui) Choose this installation if you only want the GUI and/or keep things as simple as possible.
* [**Pip installer:**](#pip) Choose this installation if you want to use AlphaViz as a Python package in an existing Python 3.8 environment (e.g. a Jupyter notebook). If needed, the GUI and CLI can be installed with pip as well.
* [**Developer installer:**](#developer) Choose this installation if you are familiar with CLI tools, [conda](https://docs.conda.io/en/latest/) and Python. This installation allows access to all available features of AlphaViz and even allows to modify its source code directly. Generally, the developer version of AlphaViz outperforms the precompiled versions which makes this the installation of choice for high-throughput experiments.

### One-click GUI

The GUI of AlphaViz is a completely stand-alone tool that requires no knowledge of Python or CLI tools. Click on one of the links below to download the latest release for:

* [**Windows**](https://github.com/MannLabs/alphaviz/releases/latest/download/alphaviz_gui_installer_windows.exe)
* [**macOS**](https://github.com/MannLabs/alphaviz/releases/latest/download/alphaviz_gui_installer_macos.pkg)
* [**Linux**](https://github.com/MannLabs/alphaviz/releases/latest/download/alphaviz_gui_installer_linux.deb)

Older releases remain available on the [release page](https://github.com/MannLabs/alphaviz/releases), but no backwards compatibility is guaranteed.

***IMPORTANT: Please refer to the [GUI manual](alphaviz/docs/alphaviz_tutorial.pdf) for detailed instructions on the installation, troubleshooting and usage of the stand-alone AlphaViz GUI.***

### Pip

AlphaViz can be installed in an existing Python 3.8 environment with a single `bash` command. *This `bash` command can also be run directly from within a Jupyter notebook by prepending it with a `!`*:

```bash
pip install alphaviz
```

Installing AlphaViz like this avoids conflicts when integrating it in other tools, as this does not enforce strict versioning of dependancies. However, if new versions of dependancies are released, they are not guaranteed to be fully compatible with AlphaViz. While this should only occur in rare cases where dependencies are not backwards compatible, you can always force AlphaViz to use dependancy versions which are known to be compatible with:

```bash
pip install "alphaviz[gui-stable]"
```

NOTE: You might need to run `pip install pip==21.0` before installing alphaviz like this. Also note the double quotes `"`.

For those who are really adventurous, it is also possible to directly install any branch (e.g. `@development`) with any extras (e.g. `#egg=alphaviz[stable,development-stable]`) from GitHub with e.g.


```bash
pip install "git+https://github.com/MannLabs/alphaviz.git@development#egg=alphaviz[stable,development-stable]"
```

### Developer

AlphaViz can also be installed in editable (i.e. developer) mode with a few `bash` commands. This allows to fully customize the software and even modify the source code to your specific needs. When an editable Python package is installed, its source code is stored in a transparent location of your choice. While optional, it is advised to first (create and) navigate to e.g. a general software folder:

```bash
mkdir ~/folder/where/to/install/software
cd ~/folder/where/to/install/software
```

***The following commands assume you do not perform any additional `cd` commands anymore***.

Next, download the AlphaViz repository from GitHub either directly or with a `git` command. This creates a new AlphaViz subfolder in your current directory.

```bash
git clone https://github.com/MannLabs/alphaviz.git
```

For any Python package, it is highly recommended to use a separate [conda virtual environment](https://docs.conda.io/en/latest/), as otherwise *dependancy conflicts can occur with already existing packages*.

```bash
conda create --name alphaviz python=3.8 -y
conda activate alphaviz
```

Finally, AlphaViz and all its [dependancies](requirements) need to be installed. To take advantage of all features and allow development (with the `-e` flag), this is best done by also installing the [development dependencies](requirements/requirements_development.txt) and/or the [gui dependencies](requirements/requirements_gui.txt) instead of only the [core dependencies](requirements/requirements.txt):

```bash
pip install -e "./alphaviz[gui,development]"
```

***By using the editable flag `-e`, all modifications to the [AlphaViz source code folder](https://github.com/MannLabs/alphaviz/tree/master/alphaviz) are directly reflected when running AlphaViz. Note that the AlphaViz folder cannot be moved and/or renamed if an editable version is installed.***

---
## Usage

There are two ways to use AlphaViz:

* [**GUI**](#gui)
* [**Python**](#python-and-jupyter-notebooks)

NOTE: The first time you use a fresh installation of AlphaViz, it is often quite slow because some functions might still need compilation on your local operating system and architecture. Subsequent use should be a lot faster.

### GUI

If the GUI was not installed through a one-click GUI installer, it can be activate with the following `bash` command:

```bash
alphaviz gui
```

Note that this needs to be prepended with a `!` when you want to run this from within a Jupyter notebook. When the command is run directly from the command-line, make sure you use the right environment (activate it with e.g. `conda activate alphaviz` or set an alias to the binary executable (can be obtained with `where alphaviz` or `which alphaviz`)).

### Python and Jupyter notebooks

AlphaViz can be imported as a Python package into any Python script or notebook with the command `import alphaviz`.

An ‘nbs’ folder in the GitHub repository contains several Jupyter Notebooks as tutorials regarding using AlphaViz as a Python package for all available pipelines: for DDA data analyzed with MaxQuant, for DIA data analyzed with DIA-NN, and for the targeted mode.

---
## Troubleshooting

In case of issues, check out the following:

* [Issues](https://github.com/MannLabs/alphaviz/issues): Try a few different search terms to find out if a similar problem has been encountered before
* [Discussions](https://github.com/MannLabs/alphaviz/discussions): Check if your problem or feature requests has been discussed before.

---
## Citations

A manuscript is currently in preparation.

---
## How to contribute

If you like this software, you can give us a [star](https://github.com/MannLabs/alphaviz/stargazers) to boost our visibility! All direct contributions are also welcome. Feel free to post a new [issue](https://github.com/MannLabs/alphaviz/issues) or clone the repository and create a [pull request](https://github.com/MannLabs/alphaviz/pulls) with a new branch. For an even more interactive participation, check out the [discussions](https://github.com/MannLabs/alphaviz/discussions) and the [the Contributors License Agreement](misc/CLA.md).

---
## Changelog

See the [HISTORY.md](HISTORY.md) for a full overview of the changes made in each version.
