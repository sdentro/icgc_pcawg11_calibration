# icgc_pcawg11_calibration
ICGC PCAWG-11 Calibration code for subclonal architecture

## Run the code
		Rscript code/1_purity_ploidy.R

		Rscript code/2_mutation_assignments.R ${sample} \
			data/morris/mutation_assignment/mutation_assignment.all.${sample}.csv \
			data/vanloo_wedge/2_clustering/${sample}/${sample}_cluster_membership.txt \
			data/peifer/Mutation_Clustering/*${sample_short}_cluster_assignments.txt \
			data/peifer/Mutation_Clustering/*${sample_short}_mclusters.txt \
			data/sahinalp/citup_pilot63/samples_v1/${sample}_cluster_membership.txt

## Dependencies - R packages
		Matrix
		RColorBrewer


## Expected directory structure

		|-- 1_purity_ploidy
		|-- 2_mutation_assignments
		|   |-- figures
		|   |-- similarities
		|   `-- tables
		|-- code
		|-- data
		|   |-- broad
		|   |-- dirichlet_input
		|   |-- morris
		|   |   |-- mutation_assignment
		|   |   `-- mutation_relations
		|   |-- peifer
		|   |   |-- Copy_Number
		|   |   |-- Mutation_Clustering
		|   |   `-- Purity_Ploidy
		|   |-- sahinalp
		|   |   `-- citup_pilot63
		|   |       `-- samples_v1
		|   `-- vanloo_wedge
		|       |-- 1_purity_ploidy
		|       |-- 2_clustering
