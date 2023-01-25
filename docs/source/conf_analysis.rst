.. _CA:

Conformation Analysis
=====================

There are many open-source tools available that enable researchers to perform conformation analyses on contact matrices. We aim to hightlight the tools and commands we use here at Cantata Bio call features. This is not a comprenisve list, nor do we own or manage any of the tools listed below. Please refer to their source document pages for more information or for help in trouble shooting the use of these tools. 
For each analysis performed below we provide, a link to the tool repo, the input file, the example command, and an example output for you to check your results against.

Example input files (.mcool and .hic files)
+++++++++++++++++++++++++++++++++++++++++++

  - Example :download:`mcool <https://s3.amazonaws.com/dovetail.pub/MicroC/mcool/MicroC_800M.mcool>` file
  - Example :download:`hic <https://s3.amazonaws.com/dovetail.pub/MicroC/hic/MicroC_800M.hic>` file


A/B Compartments
++++++++++++++++

**Recommended Software**
  
  - Fanc-C https://fan-c.readthedocs.io/en/latest/

**Example Command**

.. code-block:: console

   	fanc compartments -f -v MicroC_800M_eigen_64kb.bed -d MicroC_800M_64kb.bed -g hg38.fa MicroC_800M.mcool@64000 MicroC_800M_64kb.ab


**Example Output(s)**

  - Example :download:`AB compartments <https://s3.amazonaws.com/dovetail.pub/MicroC/AB/MicroC_800M_64kb.ab>` file
  - Example :download:`AB compartments bed <https://s3.amazonaws.com/dovetail.pub/MicroC/AB/MicroC_800M_64kb.bed>` file
  - Example :download:`AB Eigenvector <https://s3.amazonaws.com/dovetail.pub/MicroC/AB/MicroC_800M_eigen_64kb.bed>` file


Topologically Associated Domains
++++++++++++++++++++++++++++++++

**Recommended Software**

  - Juicer Arrowhead https://github.com/aidenlab/juicer 

**Example Command**

.. code-block:: console

   java -jar -Xmx48000m  -Djava.awt.headless=true -jar juicer_tools.jar arrowhead --threads 16 -k KR -m 2000 -r 10000 MicroC_800M.hic TAD_calls
   java -jar -Xmx48000m  -Djava.awt.headless=true -jar juicer_tools.jar arrowhead --threads 16 -k KR -m 2000 -r 5000 MicroC_800M.hic TAD_calls


**Example Output(s)**

  - Example :download:`10kb TADs <https://s3.amazonaws.com/dovetail.pub/MicroC/TADs/MicroC_800M_10000_blocks.bedpe>` file
  - Example :download:`5kb TADs <https://s3.amazonaws.com/dovetail.pub/MicroC/TADs/MicroC_800M_5000_blocks.bedpe>` file


Chromatin Loops
+++++++++++++++

**Recommended Software**

  - Mustache https://github.com/ay-lab/mustache 

**Example Command**

.. code-block:: console

   mustache -p 48 -f MicroC_800M.mcool -r 16000 -o MicroC_800M_16000kb_loops.tsv
   mustache -p 48 -f MicroC_800M.mcool -r 4000 -o MicroC_800M_4000kb_loops.tsv

**Example Output(s)**

  - Example :download:`16kb Loops <https://s3.amazonaws.com/dovetail.pub/MicroC/loops/MicroC_800M_16000kb_loops.tsv>` file
  - Example :download:`4kb Loops <https://s3.amazonaws.com/dovetail.pub/MicroC/loops/MicroC_800M_4000kb_loops.tsv>` file
