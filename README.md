# *C*onnected *L*inkage and *A*lignment *R*edefinition of *C*OGs (*CLARC*)

A tool that uses sequence identity, within-population linkage patterns and functional annotations to identify and reduce the over-splitting of genes into multiple clusters of orthologous genes (COGs) in a pangenome analysis. 

## Table of Contents
1. [Introduction](#introduction)
2. [Installation](#installation)
3. [Quickstart](#quickstart)
   - [Basic usage](#basic-usage)
   - [Test data](#test-data)
5. [Full usage](#full-usage)
6. [Output](#output)
7. [Raising issues](#raising-issues)
8. [Citing](#citing)

## Introduction

🚧🚧 Description under construction, come back later 🚧🚧

## Installation

### Requirements

For now, CLARC can only be installed from source by copying this github repository. Also, all dependencies are installed using a conda environment. So, the installation assumes that you have git and conda installed. 

For quick and simple instructions of how to locally install conda with miniconda follow [this link](https://docs.anaconda.com/free/miniconda/#quick-command-line-install). 

Then you can simply install git with:

```bash
conda install git
```

### Install from source

Clone repository from GitHub:

```bash
git clone https://github.com/IndraGonz/CLARC.git
```

The folder will be copied into your working directory, so make sure to navigate to your desired path for the installation.

Now, navigate to the cloned folder:

```bash
cd CLARC
```

Create conda environment with the necessary dependencies:

```bash
conda env create --file envs/clarc_env.yml
```

Activate environment:

```bash
conda activate clarc_env
```

Now from the active clarc environment, in the CLARC folder, run the installation command:

```bash
python setup.py install
```

That should do it, now you should have your very own CLARC installation within the ```clarc_env``` conda environment!

To confirm that the installation worked type:

```bash
clarc --version
```

If the version appears, then the installation was successful!

## Quickstart

### Basic usage

### Test data

## Full usage

## Output

## Raising issues

If you have trouble with any of the steps shown here, please let me know! You can email me at igonzalezojeda@g.harvard.edu or submit an issue through GitHub.

## Citing

🚧🚧 Coming soon 🚧🚧

