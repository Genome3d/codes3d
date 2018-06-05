#/bin/bash
apt-get install -y python-pip bedtools libxslt1-dev libxml2-dev
pip install configparser pandas pybedtools requests biopython matplotlib
git clone https://github.com/wikipathways/wikipathways-api-client-py.git
cd wikipathways-api-client-py
pip install -e .
cd ../
#rm -rf wikipathways-api-client-py
