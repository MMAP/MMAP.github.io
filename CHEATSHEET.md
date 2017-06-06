---
title: "MMAP Cheat Sheet"
description: "MMAP is a comprehensive mixed model program for analysis of pedigree and population data."
---
<!-- This file provides the CONTENT for the MMAP website CHEATSHEET page -->
<!-- javascript links are at the bottom of this file to improve page loading -->

<div class="toc-wrapper">
  <ol class="toc js-toc"></ol>
</div>

<p><a id="cheat_sheet" title="Cheat Sheet Intro" class="toc-item"></a></p>

# MMAP by Jeff O'Connell

[**HOME**](https://MMAP.github.io) &nbsp; &nbsp; [**Latest Release**](https://github.com/MMAP/MMAP-releases-issues-Q-and-A/releases/latest){:target="_blank"} &nbsp; &nbsp; [**Issues and Q&A**](https://github.com/MMAP/MMAP-releases-issues-Q-and-A/issues){:target="_blank"} &nbsp; &nbsp; **MMAP Cheat Sheet**

## MMAP Cheat Sheet compiled by Jim Perry

This Cheat Sheet is designed to facilitate creation of shell scripts for running MMAP.

Example soft links for the MMAP program plus input and output files are shown below.
Example commands and options using these soft links are given in subsequent sections.
Thus, you can copy the softlinks into a shell script,
adjust the link definitions, and then copy/paste the desired commands.

### Soft Links

```
mmap=/data/datasets/mmap/mmap			# path to mmap program
snps=marker_set.csv				# single column with list of markers (use SNPNAME's or RSNUM's)
subjects=subject_set.csv			# single column with list of subject id's (EGO numbers)

genoSxMbin=genotypefilename.SxM.bin		# mmap binary genotype file with Subject x Marker configuration
genoMxSbin=genotypefilename.MxS.bin		# mmap binary genotype file with Marker x Subject configuration
genoMxScsv=genotypefilename.MxS.csv		# mmap csv genotype file with Marker x Subject configuration

freqbin=genotypefilename.FREQ.mmap.bin  	# mmap binary allele frequency file
freqcsv=genotypefilename.FREQ.mmap.csv  	# mmap csv allele frequency file

output_root=myfile				# filename root for plink-like file formats (.ped, .map, .info)
``` 

<p><a id="subset_genotypes" title="Subset genotypes" class="toc-item"></a></p>

### Subset a genotype file

Creating a marker and/or subject reduced binary genotype file. The outfile is the new binary genotype file.  This command may be applied to MxS binary files (resulting in an MxS output file) or to SxM formats (resulting in SxM output file).  
`$mmap --write_reduced_genotype_binary --binary_input_filename $input_genobin --binary_output_filename $output_genobin [marker and subject options]`

#### subsetting options
```
--autosome			# Extract the autosomal SNPs
--chromosome <numbers>		# Extract SNPs on chromosomes in <numbers>  e.g.  --chromosome 1 9 22 X Y XY MT

--genomic_region <chr> <start bp> <stop bp>	#Extract SNPs in the genomic region(s) specified:  1 120000 129999

--marker_set <file>		# Extract markers in <file>
--subject_set <file>		# Extract subjects in <file>

--include_duplicate_markers`	# Use this option with `--write_reduced_genotype_binary` and `--marker_by_subject_mmap2csv`
				  to insure you get all desired markers when there are duplicate markers with SAME SNPNAME
				  are in the genotype file.
```  
**subject_set files:**  Use one column, no header, with a list of subject id's.   
**marker_set files:**  Use one column, no header, with a list of EITHER SNPNAMEs or RSNUMs.  
 &nbsp; &nbsp; ( The file type for one-column files can be .txt or .csv )

When marker_set files are only one column, then they don't need a header.
The one column can be SNPNAMEs or RSNUMs or even SNPNAMEs for some lines and RSNUMs for other lines.
MMAP will read the line, look for a match with SNPNAME in the genotype file and, if not found,
will look for a match with RSNUM in the genotype file.  MMAP seems to do everything possible to find the SNP.
Thus, if the genotype file is populated with good RSNUM data, then you can use rsNumbers in the marker_set file.
Genotype files always have SNPNAMEs, so SNPNAMEs are the best bet.

#### marker_set files can also be 2 columns (for non-sparse genotype files).
* Place SNPNAME in 1st column, and RSNUM in second column.
* If MMAP does not find the SNPNAME in the genotype file, it will look for the RSNUM (from second column).
* If the 2-column file is tab-delimited, you don't need a header (it is optional, works either way).
* If the 2-column file is comma-delimited, USE A HEADER that says: SNPNAME,RSNUM.
* If you fail to use a header with a comma-delimited file, the first row be ignored.

---

<p><a id="combine_genotypes" title="Combine genotypes" class="toc-item"></a></p>

### Combine genotype files

If chromosome-specific genotype files are created, they can then be combined into a single MMAP binary.  
`$mmap --combine_binary_genotype_files <file1> <file2> … <fileN>  --binary_output_filename <file>`

In the above command, the input and output files are MMAP marker-by-subject binary genotype files.  
RECOMMENDATION: To reduce the size of the combined binary genotype file, once the chromosome specific files are created, run the allele frequency option. The output will contain the minor allele frequency and imputation quality score, which can be used to extract a marker set based on minor allele frequency and/or imputation quality threshold. This marker set can then be used to create reduced binary genotype files before combining to the full file.

---

### Transpose a genotype file
To transpose the file from marker-by-subject (MxS) to subject-by-marker (SxM) or subject-by-marker (SxM) to marker-by-subject (MxS).   If you transpose twice then you will get the original file.  
MxS format is useful for GWA, computation of LD matrices and allele frequency calculations.  
SxM format is useful for computation of genetic covariance matrices and haplotype analysis.  

`$mmap --transpose_binary_genotype_file --binary_input_filename $genoMxSbin --binary_output_filename $genoSxMbin`

---

<p><a id="extract_genotypes" title="Extract genotypes" class="toc-item"></a></p>

### Genotype extract as csv file

The infile must be MxS format. The outfile has the results. Marker and Subjects options are valid with this command.  
`$mmap --marker_by_subject_mmap2csv --binary_input_filename $genoMxSbin --csv_output_filename $genoMxScsv`

`--include_duplicate_markers`  Use this option with `--write_reduced_genotype_binary` and `--marker_by_subject_mmap2csv` to insure you get all desired markers when there are duplicate markers with SAME SNPNAME are in the genotype file.

For spare input file use the folloing syntax. Marker_set and (soon subject_set) options are valid.  
`$mmap --sparse2csv --binary_input_filename $genoMxSbin --csv_output_filename $genoMxScsv`

---

<p><a id="extract_freq" title="Extract allele freq" class="toc-item"></a></p>

### Allele frequency extract

First create binary allele frequency file, then convert to csv.  
The allele frequency for each SNP will be in the allele frequency files.

`$mmap --write_binary_allele_frequency_file --binary_input_filename $genoMxSbin --binary_output_filename $freqbin`  
`$mmap --allele_freq_binary2csv --binary_input_filename $freqbin --csv_output_filename $freqcsv`

Alternative allele frequency extract from sparse binary:  
`$mmap --mmap_sparse2csv_allele_frequency --binary_input_filename <sparse>  --csv_output_filename <csv> `

Versions prior to 2017_03_06 used:  
`$mmap --mmap_vcf2csv_allele_frequency (even though input was a sparse file)`

---

<p><a id="haploview_files" title="Haploview files" class="toc-item"></a></p>

### Haploview extract

Creates .ped and .info file for input to Haploview program.

`$mmap --subject_by_marker_mmap2haploview --binary_input_filename $genoSxMbin --marker_set $snps --haploview_output_filename $output_root`

OPTION:  `--use_snpname` will use the SNPNAME for the variant id in the haploview dataset (.info file). Default is to use RSNUM for the variant id.

---

<p><a id="analysis_standard" title="Analysis example" class="toc-item"></a></p>

### Analysis with standard (non-imputed) genotypes ###

```
mmap=/data/datasets/mmap/mmap
ped=/data/datasets/mmap/amish.pedigree.csv 
kinbin=/data/datasets/mmap/amish.relationship.bin 

genobin=/data/datasets/markers/exmChipAdj3836.MxS.bin
pheno=/data/jperry/phenotypes/some_phenotypes.csv

covariates="sex age Exm654050 Exm654042"
suffix=exmChip
stdoutfile=z.mmap.stdout

# Clear the stdout file
cat /dev/null > $stdoutfile

# Run the regression
for trait in `cat TRAIT_LIST.txt`; 	# TRAIT_LIST.txt is a one-column file of trait names
					# Note the BACKTICK characters (they are NOT single quotes!)
do 
echo "$trait"
$mmap --ped $ped --model add --read_binary_covariance_file $kinbin \
  --phenotype_filename $pheno  --binary_genotype_filename $genobin \
  --covariates $covariates  --trait $trait  --file_suffix  $suffix \
  --binary_covariate_filename $genobin \
  --marker_set $snps \
  --subject_set $subjects \
  --min_minor_allele_frequency 0.02 \  <== ONLY for standard genotype files (see "imputation files" below)
  >> $stdoutfile
done
exit
```

---

<p><a id="analysis_example_2" title="Analysis - imputed files" class="toc-item"></a></p>

### Example for IMPUTED genotype files

$mmap --ped $ped --model add --read_binary_covariance_file $kinbin \
  --phenotype_filename $pheno  --binary_genotype_filename $genobin \
  --covariates $covariates  --trait $trait  --file_suffix  $suffix \
  --marker_set $snps \
  --subject_set $subjects \
  --output_marker_attribute INFO \  <== Add INFO column to the ...mle.pval.csv file
  --min_imputation_quality 0.3		<== minimum value for INFO 
  --chromosome 1 \					<== no leading zero on single-digit chromosome numbers
  --min_dosage 0.02 \               <== for imputation genotype files (equiv. to --min_minor_allele_frequency)
  --all_output \
  >> $stdoutfile
  
* NOTE on `--min_dosage` and `--min_minor_allele_frequency`   The "ideal" is to apply the limit based on the frequency for the subject population actually included in a model run.  However, this may use the frequency for all subjects in the genotype file.  Jeff is checking on this issue.
  
* OTHER options for Imputation analysis (these are not usually necessary)
   `--max_h2 0.98                 <== sets limit on h2 (if h2 goes to 1.0 we will get funny pValues)`
    `--snp_block_size 5000 \	  <== may not be needed (check with Jeff)`
    
 ---

<p><a id="additional_optionss" title="Additional options" class="toc-item"></a></p>

### Additional Options

If some covariates are not in your phenotype file, use one or more additional covariate files.  
MMAP will look first in the phenotype_filename, then  
 &nbsp; &nbsp; &nbsp; &nbsp; in the covariate_filename list (one or more files) and then  
 &nbsp; &nbsp; &nbsp; &nbsp; in the binary_covariate_filename list (one or more files)  
If you want to use SNPs as covariates, add the SNPname in the `--covariates` list and then  
 &nbsp; &nbsp; &nbsp; &nbsp; include the option `--binary_covariate_filename <list of one or more genotype binary files>`  
 &nbsp; &nbsp; &nbsp; &nbsp; A binary covariate file can be the same file used for `--binary_genotype_filename`  
Examples:  
```
covarFile1=/data/jperry/phenotypes/more_covariates.csv  
covarFile2=/data/jperry/phenotypes/other_covariates.csv  

--covariate_filename $covarFile1 $covarFile2
--binary_covariate_filename $genobin
```
The subject ids are expected to be in the first column of the phenotype_filename and covariate_filename files
If this condition is met, then the column header for the subject ids does not need to be specified.
If the subject ids are not in the first column, you may specify the column header as shown below
When specified, the column header must be used in phenotype_filename AND IN ALL covariate_filename files.

`--phenotype_id EGO`   # column header for subject ids in phenotype_filename and covariate_filename files.
  
---

<p><a id="gxe" title="GxE interactions" class="toc-item"></a></p>

### Gene x Environment Interactions

Gene x Environment interactions (actually, SNP x covariate interaction) are available with MMAP.
This option is valid ONLY if you are using a `--binary-genotype-filename` such that MMAP is
looping over a list of SNPs (the SNPs in the genotype file subsetted by the optional marker_set)
This option creates an additional covariate which is: SNP*covariate (SNP*BMI in example) where
SNP is the SNP from the list of SNPs being looped over.

The covariate can be any item in the phenotype file or in a covariate file or a SNP from a binary covariate file.
This item does NOT have to be in the covariate list identified with `--covariates` (but typically would be)

```
--gxe_interaction TREATMENT`  (gives additional covariate: SNP*TREATMENT)
```

NOTE: There can be only one GxE term in the model. Typically, you would have TREATMENT in the list of
covariates in addition to the --gxe_interaction term

---

<p><a id="covariate_interactions" title="Covariate interactions" class="toc-item"></a></p>

### Covariate x Covariate Interactions

You may also specify covariate x convariate interactions with the `--interactions` option.
The items do NOT have to be in the covariate list identified with `--covariates` (but typically would be)
The covariate can be any item in the phenotype file or in a covariate file or a SNP from a binary covariate file.
```
--interactions age*sex age*sex*BMI	gives 2 additional covariates: age*sex  age*sex*BMI
--interactions age*rs123456     	gives additional covariate: age*rs123456 where rs123456 is a specific SNP
					rs123456 could be in the pheno file or a covarFile or it can be in a
					binary_covariate_filename which might be the same as the
					binary_genotype_filename or could be a different binary file
```  
If you wanted to do every possible combination of 4 covariates (singles, doubles, triples, quadruples)
specify them as shown below.  If you didn't need the singles, you could leave out the `--covariates` option.
```
--covariates age sex BMI rs123  
--interactions  age*sex age*BMI age*rs123 sex*BMI sex*rs123 BMI*rs123  age*sex*BMI age*BMI*rs123 sex*BMI*rs123  age*sex*BMI*rs123
```

---

<p><a id="exclusions" title="Data exclusion" class="toc-item"></a></p>

### Exclusions based on phenotype/covariate values 

There is general syntax to exclude data from the analysis based on fields in the phenotype/covariate file.
The fields need not be used as covariates in the model itself.

examples:  
- exclude male subjects (where sex value = 1) to analyze females only. Obviously, SEX is not a covariate in the model. In the example below "SEX" must be replace by the actual header (e.g. "Gender", "sex", "SEX_abc")  
`--exclude_list 1EQ_SEX`
- AGE is a covariate in the model and we want to exclude subjects with AGE < lower_age_limit or AGE > upper_age_limit.  
`--exclude_list AGE_LT${lower_age_limit} ${upper_age_limit}LT_AGE`   
`--exclude_list AGE_LT$50 70LT_AGE    # exclude where AGE < 50 or 70 < AGE`

---

<p><a id="special_options" title="Special Options" class="toc-item"></a></p>

### Special case options:

`--polygenic_adjusted_residuals   <== Internally substitutes the "error residuals" for the trait and appends "_ADJ" to the trait name.`
`--x_male_coding_01    # treats chrX coding for males as 0 and 1  (default is 0 and 2)`

---

<p><a id="building_genotype_files" title="Building genotype files" class="toc-item"></a></p>

### Building binary genotype files

#### Create a "sparse" formated binary from vcf AND force SNPNAME to be chr:pos:ref:alt
`--vcf2mmap_binary_genotype_file --use_chr_pos_alt_ref --vcf_input_filename <vcf> --binary_output_filename <sparse>`

#### Create the MxS.bin genotype file from a .csv file

Use `--num_skip_fields 7` if we do NOT add extra columns for TYPE, GENE, AACHG  
`$mmap  --write_binary_genotype_file --csv_input_filename $genoMxScsv --binary_output_filename $genoMxSbin --num_skip_fields 7`

Use `--num_skip_fields 10` if we DO want to add 3 extra columns for TYPE, GENE, AACHG  
`$mmap  --write_binary_genotype_file --csv_input_filename $genoMxScsv --binary_output_filename $genoMxSbin --num_skip_fields 10 --additional_marker_attributes GENE C AACHG C TYPE C`

If the genotype file has the above 3 additional attributes, you must use the option
shown below (when you run mmap analysis) to have them included in the mmap output files.  
`--output_marker_attribute TYPE GENE AACHG  --snp_block_size 1`

#### Create the SxM.bin file from the MxS.bin file:  
`$mmap --transpose_binary_genotype_file --binary_input_filename $genoMxSbin --binary_output_filename $genoSxMbin`

---

<p><a id="mmap_formats" title="MMAP binary formats" class="toc-item"></a></p>

### MMAP binary formats 

- **dense** is the original MMAP binary format with file type ...bin

- **bit** format (...bit.bin) is 1/4th the size of a "dense" binary (...bin) when a binary_input_filename is required, both "bit" and "dense" formats can be used and MMAP determines the binary_type (user does not need to specify).

- **sparse** format (...sparse.bin) is the most highly compressed format.

To convert from "dense" to "bit":  
`$mmap --binary_genotype_file_dense2bit --binary_input_filename <dense> --binary_output_filename <bit>`

To convert from "bit" to "dense":  
`$mmap --binary_genotype_file_bit2dense --binary_input_filename <bit> --binary_output_filename <dense>`

To convert from "sparse" to "dense":  
`$mmap --binary_genotype_file_sparse2dense --binary_input_filename <sparse> --binary_output_filename <dense>`

To convert from "sparse" to "bit":  (NOTE: command say "dense", but we add --use_bit_coding )  
`$mmap  --binary_genotype_file_sparse2dense --use_bit_coding --binary_input_filename <sparse> --binary_output_filename <bit>`

---

<p><a id="plink_to_mmap" title="Plink to MMAP" class="toc-item"></a></p>

### Plink to MMAP binary

Convert from Plink binary to MMAP binary (assuming MxS in the Plink file):  
`$mmap --plink_bfile2mmap --swap_A1_A2 --plink_bfile $plinkBinaryFormat --binary_output_prefix $mmapFormat.MxS`  
By default, MMAP will set Plink A1 to the NON_CODED_ALLELE and A2 to the EFFECT_ALLELE  
However, Plink (by default) sets A1 to the allele with the lower allele frequency.  
If you want the resulting MMAP file to have the Plink A1 allele in MMAP's EFFECT_ALLELE,  
use: `--swap_A1_A2` which will use  
 - Plink A1 for the MMAP binary "EFFECT_ALLELE" and  
 - Plink A2 for the  "NON_CODED_ALLELE"

MMAP imports Plink binary format files into an SxM or MxS genotype binary file, depending on the Plink format, which is automatically detected.  
`$mmap --plink_bfile2mmap -–plink_bfile <prefix> --binary_output_prefix <mmap prefix>`  
Converts files \<prefix\>.bim, \<prefix\>.bed, \<prefix\>.fam into binary genotype file \<mmap prefix\>.bin and MMAP pedigree \<mmap prefix\>.ped.csv extracted from the \<prefix\>.fam.

---

<p><a id="plink_datasets" title="MMAP to Plink" class="toc-item"></a></p>

### MMAP to Plink datasets (.ped & .map)

`$mmap --subject_by_marker_mmap2plink --binary_input_filename $genoSxMbin --plink_output_prefix $output_root`  
OPTION: `--use_snpname`  will use the SNPNAME for the variant id in the plink dataset (.map file).  Default is to use RSNUM for the variant id.

`$mmap --subject_by_marker_mmap2plink --binary_input_filename <SxM binary genotype file> --plink_output_prefix <prefix>`  
Creates \<prefix\>.map and \<prefix\>.ped which can then be converted into Plink binary format with Plink commands.
Currently no support of export directly into Plink binary format.

`$mmap --marker_by_subject_mmap2tped --binary_input_filename <MxS binary genotype file> --plink_output_prefix <prefix>`  
Creates \<prefix\>.fam, \<prefix\>.bim and \<prefix\>.tped

`$mmap --marker_by_subject_mmap2plink_dosage --binary_input_filename <MxS binary genotype file> --plink_output_filename <prefix>`  
NOTE: output option may change in future to be `--plink_output_prefix`  
Creates \<prefix\>.fam, \<prefix\>.map and \<prefix\>.ped  
 - The .ped file is actually a "dosage" file  ( file type may change to ".dose" in the future )  
 - Use this option to convert an MMAP inputation binary to a plink "dosage" dataset  
NOTE: the dosage file created by this MMAP option is "format=1"  
See http://pngu.mgh.harvard.edu/~purcell/plink/dosage.shtml  
Example plink command to read files created by MMAP and create another (equivalent) dosage file.
`plink --dosage $prefix.dose format=1 --fam $prefix.fam --map $prefix.map --write-dosage --out newOut`

---

<p><a id="grm_eigenvectors" title="GRM & Eigenvectors" class="toc-item"></a></p>

### Genomic Relationship Matrices and Eigenvector files 

Eigenvector files - (to be updated soon) The eigen.bin file only depends on the subjects in the file. It is independent of trait and covariates. Create them by trait as each trait typicallly has a different number of subjects. You can add any covariate to the model as long as there is no missing data for the subjects used to create the eigenvector file.

Genomic Relationship Matrix - (to be added soon)

---

<p><a id="pheno_transforms" title="Transforming Phenotypes" class="toc-item"></a></p>

### Transformations of Phenotypes 

$mmap --ped $ped --read_binary_covariance_file $kinbin --phenotype_filename $pheno --trait $trait --all_output --covariates $covariates   --file_suffix $suff --subject_set $subject_set --transform_analysis_phenotype --binary_genotype_filename $genomxs

Takes full set of options. Exits once the file is created. NO adjustment for pedigree, which is standard. File contains all data in the model, similar to poly.model file

The pedigree, covariates and genotypes are not used in the calculation except to create a squared off data set given all the data and also to have a ready made phenotype file for analysis. 

Squared off means the data is not missing at any inputs in the model, be it pedigree, covariates and/or genotypes. When you create the transformed values, you specify the model (with pedigree, covariates, genotypes) that you plan to use the transformed values in.  For example, if you use the full phenotype file of 4700 subjects and no genotype file, then you will get a transformed file with 4700 subjects.  With a genotype file of 1100 subjects, then the squared off data set would have 1100 subjects.  The transformed values for a given subject will likely be different if the transform is done with 4700 subjects vs 1100 subjects.  If you were to take transformed data based on 4700 subjects and simply pull out the data for 1100 subjects, the histogram for those 1100 subjects might not look very normalized.

If the analysis plan calls for inverse normal of the covariate adjusted residuals, you would need 3 steps. Run the linear regression to generate the residuals, then transform the residuals with no covariates, then run mixed model.

---

<p align="center">MMAP: Mixed Model Analysis for Pedigrees and Populations - Copyright © 2017</p>
&nbsp;

<!-- And now for the javascript... -->
  <script type="text/javascript" src="https://ajax.googleapis.com/ajax/libs/jquery/1.7.1/jquery.min.js"></script>
  <script type="text/javascript" src="/assets/js_custom/application.js"></script>
