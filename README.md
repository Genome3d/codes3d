# CoDeS3D: Contextualising Developmental SNPs in Three Dimensions

## Installation

Clone the repository using the following command:

```
git clone https://github.com/alcamerone/codes3d.git
```

You will then need to install the CoDeS3D dependencies by running `setup.sh` (you may need to run this using `sudo`). CoDeS3D depends on the following:
- Ubuntu 14.04+ (untested on other systems)
- Python 2.7.9+ or Python 3+
- apt packages:
  - bedtools
  - libxslt1-dev
  - libxml2-dev
- pip packages:
  - configparser
  - pandas
  - pybedtools
  - requests
  - biopython
  - matplotlib
  - progressbar
  - progressbar2
  - psutil (run 'pip2 install psutil' for python2)
  - rpy2 (run 'pip2 install rpy2==2.8..6' for python2)
- The [WikiPathways Python API client](https://github.com/wikipathways/wikipathways-api-client-py)*

\* After installing the WikiPathways Python API client, you may get an ImportError when you attempt to run CoDeS3D. To prevent this, 
  go to the wikipathways-api-client-py/wikipathways_api_client/ directory and edit the first line in the __init__.py as follows:

`import wikipathways_api_client`

You will then be able to run CoDeS3D scripts!

### Write about miniconda set up

You will also need some data files with which to run your data. The simplest way to acquire these are using the `download_default_data.py` script. This will download the default data to the library directory specified by `libdir` in `docs/codes3d.conf`. Please note that most of the data files are very large, and this step is likely to take a long time, particularly if you wish to use all of the available HiC datasets. Note also that the total size of all available datasets is hundreds of gigabytes, and will require a large disk with a lot of free space.

## Basic Usage

The CoDeS3D interface is heavily inspired by the Qiime interface (J Gregory Caporaso *et al*., Nature Methods, 2010; doi:10.1038/nmeth.f.303). Running the CoDeS3D script in the codes3d directory will drop the user into the CoDeS3D shell, in which all CoDeS3D scripts are accessible from anywhere in the system, e.g.

```
/home/alcamerone/Documents/codes3d$ ./CoDeS3D
 Setting up CoDeS3D environment.
 
 Type 'exit' or press Ctrl+D at any time to leave.
/home/alcamerone/Documents/codes3d$ CoDeS3D> cd ../project
/home/alcamerone/Documents/project$ CoDeS3D> codes3d.py -h
 usage: codes3d.py [-h] -i INPUTS [INPUTS ...] [-c CONFIG]
                  [-n INCLUDE_CELL_LINES [INCLUDE_CELL_LINES ...]]
                  [-x EXCLUDE_CELL_LINES [EXCLUDE_CELL_LINES ...]]
                  [-o OUTPUT_DIR] [-l] [-s] [-p NUM_PROCESSES]
                  [-f FDR_THRESHOLD]
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
