.. _BYB:

Before you begin
================

Have a copy of the Micro-C scripts on your machine:
----------------------------------------------------

Clone this repository:

.. code-block:: console

   git clone https://github.com/dovetail-genomics/Micro-C.git


Dependencies
-------------

Make sure that the following dependencies are installed:

- `pysam <https://pysam.readthedocs.io/en/latest/>`_
- `tabulate <https://pypi.org/project/tabulate/>`_
- `bedtools <https://bedtools.readthedocs.io/en/latest/index.html>`_
- `deeptools <https://deeptools.readthedocs.io/en/develop/>`_
- `matplotlib <https://matplotlib.org/>`_
- `pandas <https://pandas.pydata.org/pandas-docs/stable/dsintro.html>`_
- `bwa <https://github.com/lh3/bwa>`_
- `pairtools <https://github.com/open2c/pairtools>`_
- `preseq <http://smithlabresearch.org/software/preseq/>`_
- `samtools <https://github.com/samtools/samtools>`_

If you are facing any issues with the installation of any of the dependencies, please contact the supporter of the package.

python3 and pip3 are required, if you don't already have them installed, you will need sudo privileges.

* Update and install python3 and pip3:

.. code-block:: console 

   sudo apt-get update
   sudo apt-get install python3 python3-pip


* To set python3 and pip3 as primary alternative:

.. code-block:: console

   sudo update-alternatives --install /usr/bin/python python /usr/bin/python3 1
   sudo update-alternatives --install /usr/bin/pip pip /usr/bin/pip3 1


If you are working on a new machine and don't have the dependencies, you can use the ``installDep.sh`` script in this repository for updating your instance and installing the dependencies and python3. This process will take approximately 10' and requires sudo privileges. The script was tested on Ubuntu 18.04 with the latest version as of 04/11/2020

If you choose to run the provided installation script you will first need to set the permission to the file:

.. code-block:: console

   chmod +x ./Micro-C/installDep.sh


And then run the installation script:

.. code-block:: console

   ./Micro-C/installDep.sh


.. admonition:: Remember!

   Once the installetion is completed, sign off and then sign back to your instance to refresh the database of applications.


Input files
-----------

For this tutorial you will need: 

* **fastq files** R1 and R2, either fastq or fastq.gz are acceptable
* **reference in a fasta file format**, e.g. hg38

If you don't already have your own input files or want to run a test on a small data set, you can download sample fastq files from the :ref:`Micro-C Data Sets section<DATASETS>`. The 2M data set is suitable for a quick testing of the instruction in this tutorial. 

.. code-block:: console

   wget https://s3.amazonaws.com/dovetail.pub/HiC/fastqs/MicroC_2M_R1.fastq
   wget https://s3.amazonaws.com/dovetail.pub/HiC/fastqs/MicroC_2M_R2.fastq

