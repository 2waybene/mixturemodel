#!/bin/sh

export QIIME_CONFIG_FP=/home/qglo4/Analysis_tools/Qiime/qiime-1.7.0-release/../qiime_config
export LAPACK=None
export PYTHONPATH=/home/qglo4/Analysis_tools/Qiime/emperor-0.0.0-repository-245e0563/lib/python2.7/site-packages:/home/qglo4/Analysis_tools/Qiime/emperor-0.0.0-repository-245e0563/lib/:/home/qglo4/Analysis_tools/Qiime/qiime-galaxy-0.0.1-repository-5b980770/lib/:/home/qglo4/Analysis_tools/Qiime/qiime-1.7.0-release/lib/python2.7/site-packages:/home/qglo4/Analysis_tools/Qiime/qiime-1.7.0-release/lib/:/home/qglo4/Analysis_tools/Qiime/matplotlib-1.1.0-release/lib/python2.7/site-packages:/home/qglo4/Analysis_tools/Qiime/pprospector-1.0.1-release/lib/python2.7/site-packages:/home/qglo4/Analysis_tools/Qiime/tax2tree-1.0-release/lib/python2.7/site-packages:/home/qglo4/Analysis_tools/Qiime/pynast-1.2-release/lib/python2.7/site-packages:/home/qglo4/Analysis_tools/Qiime/pycogent-1.5.3-release/lib/python2.7/site-packages:/home/qglo4/Analysis_tools/Qiime/biom-format-1.1.2-release/lib/python2.7/site-packages:/home/qglo4/Analysis_tools/Qiime/sphinx-1.0.4-release/lib/python2.7/site-packages:/home/qglo4/Analysis_tools/Qiime/MySQL-python-1.2.3-release/lib/python2.7/site-packages:/home/qglo4/Analysis_tools/Qiime/SQLAlchemy-0.7.1-release/lib/python2.7/site-packages:/home/qglo4/Analysis_tools/Qiime/mpi4py-1.2.2-release/lib/python2.7/site-packages:/home/qglo4/Analysis_tools/Qiime/setuptools-0.6c11-release/lib/python2.7/site-packages:/home/qglo4/Analysis_tools/Qiime/gdata-2.0.17-release/lib/python2.7/site-packages:/home/qglo4/Analysis_tools/Qiime/ipython-latest-repository-08f5c195/lib/python2.7/site-packages:/home/qglo4/Analysis_tools/Qiime/tornado-2.2.1-release/lib/python2.7/site-packages:/home/qglo4/Analysis_tools/Qiime/numpy-1.5.1-release/lib/python2.7/site-packages:/home/qglo4/Analysis_tools/Qiime/cython-0.17-release/lib/python2.7/site-packages:/home/qglo4/Analysis_tools/Qiime/pyzmq-2.1.11-release/lib/python2.7/site-packages:$PYTHONPATH
export BLAS=None
export TEST_DB=1
export BLASTMAT=/home/qglo4/Analysis_tools/Qiime/blast-2.2.22-release/data
export ATLAS=None
export PYRO_LOOKUP_FILE=/home/qglo4/Analysis_tools/Qiime/ampliconnoise-1.27-release/Data/LookUp_E123.dat
export RDP_JAR_PATH=/home/qglo4/Analysis_tools/Qiime/rdpclassifier-2.2-release/rdp_classifier-2.2.jar
export PYCOGENT=/home/qglo4/Analysis_tools/Qiime/pycogent-1.5.3-release/./
export SEQ_LOOKUP_FILE=/home/qglo4/Analysis_tools/Qiime/ampliconnoise-1.27-release/Data/Tran.dat
export PATH=/home/qglo4/Analysis_tools/Qiime/chimeraslayer-4.29.2010-release/ChimeraSlayer:/home/qglo4/Analysis_tools/Qiime/chimeraslayer-4.29.2010-release/NAST-iEr:/home/qglo4/Analysis_tools/Qiime/emperor-0.0.0-repository-245e0563/bin:/home/qglo4/Analysis_tools/Qiime/qiime-galaxy-0.0.1-repository-5b980770/scripts:/home/qglo4/Analysis_tools/Qiime/qiime-1.7.0-release/bin:/home/qglo4/Analysis_tools/Qiime/pprospector-1.0.1-release/bin:/home/qglo4/Analysis_tools/Qiime/tax2tree-1.0-release/bin:/home/qglo4/Analysis_tools/Qiime/pynast-1.2-release/bin:/home/qglo4/Analysis_tools/Qiime/biom-format-1.1.2-release/bin:/home/qglo4/Analysis_tools/Qiime/cython-0.17-release/bin:/home/qglo4/Analysis_tools/Qiime/ipython-latest-repository-08f5c195/bin:/home/qglo4/Analysis_tools/Qiime/tornado-2.2.1-release/bin:/home/qglo4/Analysis_tools/Qiime/rdpclassifier-2.2-release/.:/home/qglo4/Analysis_tools/Qiime/pyzmq-2.1.11-release/bin:/home/qglo4/Analysis_tools/Qiime/python-2.7.3-release/bin:/home/qglo4/Analysis_tools/Qiime/vienna-1.8.4-release/.:/home/qglo4/Analysis_tools/Qiime/uclust-1.2.22-release/.:/home/qglo4/Analysis_tools/Qiime/cdhit-3.1-release/.:/home/qglo4/Analysis_tools/Qiime/infernal-1.0.2-release/bin:/home/qglo4/Analysis_tools/Qiime/parsinsert-1.0.4-release/.:/home/qglo4/Analysis_tools/Qiime/muscle-3.8.31-release/.:/home/qglo4/Analysis_tools/Qiime/bwa-0.6.2-release/.:/home/qglo4/Analysis_tools/Qiime/blast-2.2.22-release/bin:/home/qglo4/Analysis_tools/Qiime/pplacer-1.1-release/.:/home/qglo4/Analysis_tools/Qiime/fasttree-2.1.3-release/.:/home/qglo4/Analysis_tools/Qiime/cdbtools-10.11.2010-release/.:/home/qglo4/Analysis_tools/Qiime/mothur-1.25.0-release/.:/home/qglo4/Analysis_tools/Qiime/cytoscape-2.7.0-release/.:/home/qglo4/Analysis_tools/Qiime/drisee-1.2-release/.:/home/qglo4/Analysis_tools/Qiime/ampliconnoise-1.27-release/Scripts:/home/qglo4/Analysis_tools/Qiime/ampliconnoise-1.27-release/bin:/home/qglo4/Analysis_tools/Qiime/raxml-7.3.0-release/.:/home/qglo4/Analysis_tools/Qiime/sourcetracker-0.9.5-release/.:/home/qglo4/Analysis_tools/Qiime/rtax-0.983-release/.:/home/qglo4/Analysis_tools/Qiime/clearcut-1.0.9-release/.:/home/qglo4/Analysis_tools/Qiime/blat-34-release/.:/usr/local/bin:/usr/bin:/bin:$PATH
export SOURCETRACKER_PATH=/home/qglo4/Analysis_tools/Qiime/sourcetracker-0.9.5-release/.
export LD_LIBRARY_PATH=/home/qglo4/Analysis_tools/Qiime/python-2.7.3-release/lib
export QIIME=/home/qglo4/Analysis_tools/Qiime/qiime-1.7.0-release/./



cd /home/qglo4/Analysis_tools/QAAP/Code/Analysis_V_2.1
perl start_qiime_pipeline.pl




