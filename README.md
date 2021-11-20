# Merge and Traverse Focus Assembler

Next Generation Sequencing (NGS) technologies are capable of producing millions to even billions of DNA fragments called reads. These reads are assembled into larger sequences called contigs by graph theoretic software tools called assemblers. The Focus Assembler is a multilayer graph model for the analysis and extraction of biologically relevant features in NGS data sets. Unlike previous models that use a single graph to model the assembly graph, Focus utilizes a multiset of graphs across a spectrum of granularity. A single graph is only capable of providing a single view of the read data set; however, a multset of graphs is able to capture both local and global information at different scales of granularity. 

## References
### General Algorithm

Warnke-Sommer, J., Ali, H. Graph mining for next generation sequencing: leveraging the assembly graph for biological insights. BMC Genomics 17, 340 (2016). https://doi.org/10.1186/s12864-016-2678-2 (https://doi.org/10.1186/s12864-016-2678-2)

### Efficiency and Succinct Data Structures

J. Warnke and H. H. Ali, "An efficient overlap graph coarsening approach for modeling short reads," 2012 IEEE International Conference on Bioinformatics and Biomedicine Workshops, 2012, pp. 704-711, doi: 10.1109/BIBMW.2012.6470223 (https://doi.org/10.1109/BIBMW.2012.6470223)

### HPC Parallel Sequence Alignment 

J. Warnke, S. Pawaskar and H. Ali, "An energy-aware bioinformatics application for assembling short reads in high performance computing systems," 2012 International Conference on High Performance Computing & Simulation (HPCS), 2012, pp. 154-160, doi: 10.1109/HPCSim.2012.6266905 (https://doi.org/10.1109/HPCSim.2012.6266905)

### HPC Distributed Hybrid Graph 

J. D. Warnke-Sommer and H. H. Ali, "Parallel NGS Assembly Using Distributed Assembly Graphs Enriched with Biological Knowledge," 2017 IEEE International Parallel and Distributed Processing Symposium Workshops (IPDPSW), 2017, pp. 273-282, doi: 10.1109/IPDPSW.2017.143 (https://doi.org/10.1109/IPDPSW.2017.143)

## Dependencies 

* openmpi (https://www.open-mpi.org/)

## Installation

``` cd focus_assembler ```

``` make clean ```

Making the serial components 

``` make all ```

Making the parallel components

``` make -f Makefile.MPI ```
	
The resulting executables will be in the bin directory.

## Usage

The Focus assembler is composed of 4 required modules and a fifth optional module. 

1. Preprocessing 

	This module trims the NGS reads and removes low quality sequence.

	```/Path/To/Dir/focus_assembler/bin/preprocess --workDir workDir --numJobs 4 --singleReads /Path/To/Dir/focus_assembler/data/test.fasta``` 

2. Alignment 

	This HPC module performs parallel alignment of the NGS reads. 

	``` mpirun /Path/To/Dir/focus_assembler/bin/parallelAlign --workDir workDir  --minIden 90 --minOverlap 50 ```

3. Coarsen Graph

	This module creates the coarsened graph spectrum.

	``` mpirun /Path/To/Dir/focus_assembler/bin/coarsen --workDir workDir --minDensity 5 --percentMerged .001 --verbose --time ```

4. UnCoarsen Graph 

	This module creates the hybrid graph spectrum.

	``` mpirun /Path/To/Dir/focus_assembler/bin/unCoarsen --workDir wordDir --time ```
 
5. Distributed Graph Spectrum (optional)

	This module implements the Kernighan-Lin algorithm to partition the hybrid graph. 

	``` mpirun /Path/To/Dir/focus_assembler/bin/distributeSpectrum --workDir wordDir --time --distributedDir distribDir --numPartitions 4```


### Additional Notes
Script examples for submiting assembly jobs to the SLURM job scheduler can be found in the scripts folder. 

A test data set, (test.fasta, test.qual), can be found in the data directory.

Run any program without arguements to display a help page with program arguments.

## License
