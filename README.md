# Block_bootstrap_for_vcf
This pipeline performs block bootstrap resampling of a VCF for downstream demographic inference (e.g. with GADMA).

Workflow:

Extract scaffold lengths from the VCF header using bcftools.
Split scaffolds into blocks based on a recombination-derived distance (R / data.table).
Generate bootstrap replicates by sampling blocks with replacement.
Reconstruct bootstrap VCFs (preserving duplicate regions) with bcftools parallely.
Compute SFS for each bootstrap using easySFS parallely.

Input:

VCF file
Popmap file
Projection values

Output:

Bootstrap VCFs
Corresponding SFS files (one per bootstrap replicate)
