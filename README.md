# MetaPathways for FOGDOG

## Installation

MetaPathways for FOGDOG requires: 
*Python 2.7 or greater 
*Scipy (python)
*python biom-format package
*[Pathway Tools](http://bioinformatics.ai.sri.com/ptools/) developed by SRI International

In addtion MetaPathways uses several third party softwares all of which are distributed with the MetaPathways code. Licenses for each third party software are within the tar ball (source code). Note that these third party softwares have all been modified for compatability with MetaPathways as is permissable under the included licenses. 

After the repositpry in cloned, navigate into MetaPathways_Python_Koonkie.3.0/executables/source In this folder you will see a Makefile. The following are options:

make untar_folder
make build_folders
make install_macosx
make install_redhat
make install_ubuntu

it is also possible to:
make clean_folders
make tar_folders
make remove_macosx
make remove_redhat
make remove_ubuntu

Please see the [MetaPathways v2.5 wiki](https://github.com/hallamlab/metapathways2/wiki) for more details on running the software.

A template [MetaPathways_DBs.zip (**Updated: October 2014**)](https://www.dropbox.com/s/ye3kpve041e0r39/MetaPathways_DBs.zip?dl=0) contains starter protein and taxonomic databases

## Citation

This version of MetaPathways is a modified version of Metapathway v2.0. It is modified as is permissable according to the license. The citation and abstract can be found below.

Niels W. Hanson, Kishori M. Konwar, Shang-Ju Wu, Steven J. Hallam. *MetaPathways v2.0: A master-worker model for environmental Pathway/Genome Database construction on grids and clouds.* Proceedings of the 2014 IEEE Conference on Computational Intelligence in Bioinformatics and Computational Biology (CIBCB 2014), Honolulu, HI, USA, May 21-24, 2014. [doi:10.1109/CIBCB.2014.6845516](http://ieeexplore.ieee.org/xpl/articleDetails.jsp?arnumber=6845516)


The development of high-throughput sequencing technologies over the past decade has generated a tidal wave of environmental sequence information from a variety of natural and human engineered ecosystems. The resulting flood of infor- mation into public databases and archived sequencing projects has exponentially expanded computational resource requirements rendering most local homology-based search methods inefficient. We recently introduced MetaPathways v1.0, a modular annotation and analysis pipeline for constructing environmental Pathway/Genome Databases (ePGDBs) from environmental sequence information capable of using the Sun Grid engine for external resource partitioning. However, a command-line interface and facile task management introduced user activation barriers with concomitant decrease in fault tolerance.

Here we present MetaPathways v2.0 incorporating a graphical user interface (GUI) and refined task management methods. The MetaPathways GUI provides an intuitive display for setup and process monitoring and supports interactive data visualization and sub-setting via a custom Knowledge Engine data structure. A master-worker model is adopted for task management allowing users to scavenge computational results from a number of worker grids in an ad hoc, asynchronous, distributed network that dramatically increases fault tolerance. This model facilitates the use of EC2 instances extending ePGDB construction to the Amazon Elastic Cloud.
