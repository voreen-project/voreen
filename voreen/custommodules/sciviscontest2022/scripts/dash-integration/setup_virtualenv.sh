#!/bin/bash

pip3 install virtualenv==16.7.12
virtualenv -p python3 voreenenv
source voreenenv/bin/activate
pip install dash
pip install scikit-learn
pip install umap-learn
pip install matplotlib
pip install pandas
#pip uninstall numpy
#pip install numpy==1.19
