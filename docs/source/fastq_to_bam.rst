.. _FTB:

From fastq to final valid pairs bam file
========================================

.. _Impatient: 

.. admonition:: fastq to final valid pairs bam file - for the impatient!

   If you just want to give it a shot and run all the alignment and filtering steps without going over all the details, we made a shorter version for you, with all the steps piped, outputting a final bam file with its index file and a dup stats file, otherwise move to the next section :ref:`fastq to final valid pairs bam file - step by step<step-by-step>`
   
   **Command:**
   
   .. code-block:: console
   
      bwa mem -5SP -T0 -t<cores> <ref.fa> <MicroC.R1.fastq.gz> <MicroC.R2.fastq.gz>| \ 
      pairtools parse --min-mapq 40 --walks-policy 5unique \ 
      --max-inter-align-gap 30 --nproc-in <cores> --nproc-out <cores> --chroms-path <ref.genome> | \ 
      pairtools sort --tmpdir=<full_path/to/tmpdir> --nproc <cores>|pairtools dedup --nproc-in <cores> \ 
      --nproc-out <cores> --mark-dups --output-stats <stats.txt>|pairtools split --nproc-in <cores> \ 
      --nproc-out <cores> --output-pairs <mapped.pairs> --output-sam -|samtools view -bS -@<cores> | \
      samtools sort -@<cores> -o <mapped.PT.bam>;samtools index <mapped.PT.bam>
   
   **Example:**
   
   .. code-block:: console
   
      bwa mem -5SP -T0 -t16 hg38.fasta MicroC_2M_R1.fastq MicroC_2M_R2.fastq| pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in 8 --nproc-out 8 --chroms-path hg38.genome | pairtools sort --tmpdir=/home/ubuntu/ebs/temp/ --nproc 16|pairtools dedup --nproc-in 8 --nproc-out 8 --mark-dups --output-stats stats.txt|pairtools split --nproc-in 8 --nproc-out 8 --output-pairs mapped.pairs --output-sam -|samtools view -bS -@16 | samtools sort -@16 -o mapped.PT.bam;samtools index mapped.PT.bam


|clock| The full command above, with 2M read pairs on an Ubuntu 18.04 machine with 16 CPUs and 64GiB was completed in less than 5 minutes. 
On the same machine type.

.. |clock| image:: /images/clock.jpg
           :scale: 5 %

.. _step-by-step:

fastq to final valid pairs bam file - step by step
--------------------------------------------------

Alignment 
+++++++++

Now that you have a genome file, index file and a reference fasta file you are all set to align your Micro-C library to the reference. Please note the specific settings that are needed to map mates independently and for optimal results with our proximity library reads.


.. csv-table::
   :file: tables/alignment.csv
   :header-rows: 1
   :widths: 25 75
   :class: tight-table


Bwa mem will output a sam file that you can either pipe or save to a path using -o option, as in the example below:

**Command:**

.. code-block:: console

   bwa mem -5SP -T0 -t<threads> <ref.fasta> <MicroC_R1.fastq> <MicroC_R2.fastq> -o <aligned.sam> 


**Example (one pair of fastq files):**

.. code-block:: console

   bwa mem -5SP -T0 -t16 hg38.fasta MicroC_2M_R1.fastq MicroC_2M_R2.fastq -o aligned.sam


**Example (multiple pairs of fastq files):**

.. code-block:: console

   bwa mem -5SP -T0 -t16 hg38.fasta <(cat file1.R1.fastq.gz file2.R1.fastq.gz file3.R1.fastq.gz) <(file1.R2.fastq.gz file2.R2.fastq.gz file3.R2.fastq.gz) -o aligned.sam



Recording valid ligation events
+++++++++++++++++++++++++++++++

We use the ``parse`` module of the ``pairtools`` pipeline to find ligation junctions in Micro-C (and other proximity ligation) libraries. When a ligation event is identified in the alignment file the pairtools pipeline will record the outer-most (5’) aligned base pair and the strand of each one of the paired reads into ``.pairsam`` file (pairsam format captures SAM entries together with the Hi-C pair information). In addition, it will also asign a pair type for each event. e.g. if both reads aligned uniquely to only one region in the genome, the type UU (Unique-Unique) will be assigned to the pair. The following steps are necessary to identify the high quality valid pairs over low quality events (e.g. due to low mapping quality):


``pairtools parse`` options:


.. csv-table::
   :file: tables/parse.csv
   :header-rows: 1
   :widths: 20 20 60
   :class: tight-table


``pairtools parse`` command example for finding ligation events:

**Command:**

.. code-block:: console

   pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in <cores>\
   --nproc-out <cores> --chroms-path <ref.genome> <aligned.sam> > <parsed.pairsam>


**Example:**

.. code-block:: console

   pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in 8 --nproc-out 8 --chroms-path hg38.genome aligned.sam >  parsed.pairsam


At the parsing step, pairs will be flipped such that regardless of read1 and read2, pairs are always recorded with first side of the pair having the lower genomic coordinates. 


Sorting the pairsam file
++++++++++++++++++++++++


The parsed pairs are then sorted using `pairtools sort`

``pairtools sort`` options:

.. csv-table::
   :file: tables/sort.csv
   :header-rows: 1
   :widths: 25 75
   :class: tight-table

**Command:**

.. code-block:: console

   pairtools sort --nproc <cores> --tmpdir=<path/to/tmpdir> <parsed.pairsam> > <sorted.pairsam>


**Example:**

.. code-block:: console

   pairtools sort --nproc 16 --tmpdir=/home/ubuntu/ebs/temp/  parsed.pairsam > sorted.pairsam


.. admonition:: Important!

   Please note that an absolute path for the temp directory is required for ``pairtools sort``, e.g. path of the structure ~/ebs/temp/ or ./temp/ will not work, instead, something of this sort is needed /home/user/ebs/temp/


.. _DUPs:

Removig PCR duplicates
++++++++++++++++++++++

``pairtools dedup`` detects molecules that could be formed via PCR duplication and tags them as “DD” pair type. These pairs should be excluded from downstream analysis. Use the pairtools dedup command with the `--output-stats` option to save the dup stats into a text file.

```pairtools dedup``` options:

.. csv-table::
   :file: tables/dedup.csv
   :header-rows: 1
   :widths: 25 75
   :class: tight-table

**Command:**

.. code-block:: console

   pairtools dedup --nproc-in <cores> --nproc-out <cores> --mark-dups --output-stats <stats.txt> \
   --output <dedup.pairsam> <sorted.pairsam>


**Example:**

.. code-block:: console

   pairtools dedup --nproc-in 8 --nproc-out 8 --mark-dups --output-stats stats.txt --output dedup.pairsam sorted.pairsam

.. _GPB:

Generate .pairs and bam files
+++++++++++++++++++++++++++++

The ``pairtools split`` command is used to split the final ``.pairsam`` into two files: ``.sam`` (or ``.bam``) and ``.pairs`` (``.pairsam`` has two extra columns containing the alignments from which the Micro-C pair was extracted, these two columns are not included in ``.pairs`` files)

``pairtools split`` options:

.. csv-table::
   :file: tables/split.csv
   :header-rows: 1
   :widths: 25 75
   :class: tight-table


**Command:**

.. code-block:: console

   pairtools split --nproc-in <cores> --nproc-out <cores> --output-pairs <mapped.pairs> \
   --output-sam <unsorted.bam> <dedup.pairsam>


**Example:**

.. code-block:: console

   pairtools split --nproc-in 8 --nproc-out 8 --output-pairs mapped.pairs --output-sam unsorted.bam dedup.pairsam

The ``.pairs`` file can be used for generating :ref:`contact matrix <GCM>`

Generating the final bam file
+++++++++++++++++++++++++++++

For downstream steps, the bam file should be sorted, using the command `samtools sort`

``samtools sort`` options:

.. csv-table::
   :file: tables/bam_sort.csv
   :header-rows: 1
   :widths: 25 75
   :class: tight-table
 

**Command:**

.. code-block:: console

   samtools sort -@<threads> -T <path/to/tmpdir/tempfile.bam>-o <mapped.PT.bam> <unsorted.bam>


**Example:**

.. code-block:: console

   samtools sort -@16 -T /home/ubuntu/ebs/temp/temp.bam -o mapped.PT.bam unsorted.bam


For future steps an index (.bai) of the bam file is also needed.
Index the bam file:

**Command:**

.. code-block:: console

   samtools index <mapped.PT.bam>


**Example:**

.. code-block:: console

   samtools index mapped.PT.bam


The ``mapped.PT.bam`` is the final bam file that will be used downstream steps.

The above steps resulted in multiple intermediate files, to simplify the process and avoid intermediate files, you can pipe the steps as in the example above (:ref:`fastq to final valid pairs bam file - for the impatient<Impatient>`)

