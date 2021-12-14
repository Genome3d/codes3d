# CoDeS3D: Contextualising Developmental SNPs in Three Dimensions

## Installation

Clone the repository using the following command:

```
https://github.com/Genome3d/codes3d.git
cd codes3d/
```

## Setup
(The following guide has been tested on Ubuntu 14.04+.)

Install required apt packages:
```
sudo apt update
sudo apt install -y libz-dev libxslt1-dev libpq-dev python3-pip python3-dev python3-virtualenv bedtools
pip3 install --upgrade pip virtualenv setuptools
```


#### Create environment with conda
Ensure that you have [conda](https://docs.conda.io/en/latest/) installed. Then create a conda environment from the root directory. In this example, the environment is created in the envs directory and has all the required dependencies. (Size of environment ~1.8GB.)
```
conda env create --prefix ./codes3d_env --file environment.yaml
```

Then activate the environment with
```
conda activate codes3d_env/
```
To deactivate,
```
conda deactivate
```


#### Data installation
You will also need to install Hi-C libraries and data necessary to calculate eQTLs. You can find detailed instructions and helper scripts in the `download_data` directory's [README](download_data/README.md) file.

## Basic Usage

The CoDeS3D interface is heavily inspired by the Qiime interface (J Gregory Caporaso *et al*., Nature Methods, 2010; doi:10.1038/nmeth.f.303). Running the CoDeS3D script in the codes3d directory will drop the user into the CoDeS3D shell, in which all CoDeS3D scripts are accessible from anywhere in the system, e.g.

```
/home/alcamerone/Documents/codes3d$ conda activate envs/codes3d_env
(codes3d_env) /home/alcamerone/Documents/codes3d$ ./CoDeS3D
 Setting up CoDeS3D environment.

 Type 'exit' or press Ctrl+D at any time to leave.
/home/alcamerone/Documents/codes3d$ CoDeS3D> cd ../project
/home/alcamerone/Documents/project$ CoDeS3D> codes3d.py -h
usage: codes3d.py [-h] [-s SNP_INPUT [SNP_INPUT ...]]                    
                 [-g GENE_INPUT [GENE_INPUT ...]] [-o OUTPUT_DIR]    
                 [--multi-test MULTI_TEST] [--pval-threshold PVAL_THRESHOLD]
                 [-f FDR_THRESHOLD] [--maf-threshold MAF_THRESHOLD]       
                 [-p NUM_PROCESSES] [--no-afc]                       
                 [--afc-bootstrap AFC_BOOTSTRAP]                    
                 [-n INCLUDE_CELL_LINES [INCLUDE_CELL_LINES ...]]        
                 [-x EXCLUDE_CELL_LINES [EXCLUDE_CELL_LINES ...]]
                 [--list-hic-libraries] [--match-tissues ...]                                                      
                 [--list-tissue-tags] [-t TISSUES [TISSUES ...]]          
                 [--eqtl-project EQTL_PROJECT] [--list-eqtl-db]            
                 [-r RESTRICTION_ENZYMES [RESTRICTION_ENZYMES ...]]    
                 [--list-enzymes] [--list-eqtl-tissues] [-c CONFIG]      
                 [--output-format OUTPUT_FORMAT] [--do-not-produce-summary]
                 [--suppress-intermediate-files] [--non-spatial]        
 [...]
```

If you want to make CoDeS3D accessible from anywhere in your system, you can add the codes3d directory to your PATH by appending this line to the bottom of your `.bash_profile`/`.bashrc` file in your home directory:

```
export PATH=$PATH:/path/to/CoDeS3D # Replace "/path/to" with the path where the CoDeS3D program is stored
```

Alternatively, you can make a link to the program in a directory which is already in your path (e.g. `/usr/local/bin`) using the following commands:

```
cd /usr/local/bin
sudo ln -s /path/to/CoDeS3D CoDeS3D # Replace "/path/to" with the path where the CoDeS3D program is stored
```

For details on the expected inputs of any CoDeS3D script, simply run the script with the `-h` argument, as in the example above.
