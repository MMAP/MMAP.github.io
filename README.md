<!-- This file provides the CONTENT for the MMAP website -->
<!-- javascript links are at the bottom of this file to improve page loading -->

# MMAP by Jeff O'Connell

[**Latest Release**](https://github.com/MMAP/MMAP-releases-issues-Q-and-A/releases/latest){:target="_blank"} &nbsp; &nbsp; [**Issues and Q&A**](https://github.com/MMAP/MMAP-releases-issues-Q-and-A/issues){:target="_blank"}

## MMAP: Mixed Model Analysis for Pedigrees and Populations

MMAP is a comprehensive mixed model program for analysis of pedigree and population data. It provides an optimized and flexible platform that incorporates a wide range of covariance structures such as additive, dominance, epistasis, maternal and imprinting using pedigree and/or genotype data and also allows users to define their own covariance structures. Likelihood calculations use multi-threaded optimized matrix libraries to handle multiple random effects. MMAP can import data from a variety of imputation programs to avoid data manipulation and IBS/IBD programs to build covariance structures.

MMAP uses a fast low-memory method to calculate additive and dominant genetic covariance structures using SNP data, which can be quite challenging for large data sets. For polygenic SNP analysis MMAP can store SNP-covariance products to reduce the complexity subsequent analyses with the same subjects to linear regression, independent of phenotype or covariates.

<p><a id="installation" title="Installation" class="toc-item"></a></p>

### Installation

MMAP is compiled with the Intel Math Kernel library for the Unix/Linux environment and uses BLAS and LAPACK libraries. To ensure compatibility only a static executable is currently available. After download, be sure to make the file executable. [Click here to download MMAP](https://github.com/MMAP/MMAP-releases-issues-Q-and-A/releases/latest){:target="_blank"}.

<p><a id="pedigree" title="Pedigree Files" class="toc-item"></a></p>

### Pedigree Files

`--ped <csv pedigree file>` identifies the comma delimited pedigree file with header containing the first 5 tokens:  
```
  PED  - Pedigree ID
  EGO  - Individual ID
  FA   - Father ID
  MO   - Mother ID
  SEX  - sex (1 = male, 2 = female)
```

The ID is alphanumeric. MMAP expects IDs to be unique across pedigrees as there is no pedigree information in the phenotype file. Missing parents are coded as 0. Parents that do not have records as individuals are ignored and single missing parents are allowed. Unknown sex is not supported. There are two additional columns (not ordered) that will be interpreted by MMAP, if present.
```
  MZTWIN - non-zero integer as group identifier of genetically identical individuals
  COHORT - non-zero integer defines group effect (see section on variance components)
```
  
The pedigree is assumed to be in ancestral order which means parents are listed before offspring. Ancestral order is used for efficient computation of the relationship matrix. MMAP currently enforces ancestral order but allows parents to be listed that do have their own records. If errors are encountered, MMAP exists and errors are output to the file pedigree.csv, which can used to correct the problems.

Additional pedigree option:

`--single_pedigree` instructs MMAP to ignore the pedigree ID and interpret the pedigree file as a single pedigree. This option is useful when creating covariance matrices that include between-pedigree values such as genetic similarity through genotype data.

Example <!--[pedigree file](files/pedigree.txt){:download="pedigree.txt"}-->[pedigree file](files/pedigree.txt){:target="_blank"} contains two pedigrees. In the first pedigree there are two sets of genetically identical subjects: MGM and PGF, and 1, 2 and 3.

<p><a id="phenotypes" title="Phenotype Files" class="toc-item"></a></p>

### Phenotype Files

Phenotypes are stored in a comma delimited file with header. The first column is assumed to be individual ID (independent of actual header token). Currently only numerical traits and covariates are supported; they are read in as real numbers, but can be coded as integers in the file. Missing data is coded as blank.

`--phenotype_filename <phenotype file>` identifies the phenotype file to be read in.  
`--phenotype_id <header token>` identifies column containing individual id. Not needed if individual id is the first column.

`--trait <trait1> <trait2> ... <traitN>` identifies the traits to be analyzed.  
 &nbsp; NOTE: Currently only a single trait is supported. Multitrait analysis is being implemented.

Example [phenotype file](files/phenotype.txt){:target="_blank"} contains data for 5 subjects from the [pedigree file](files/pedigree.txt){:target="_blank"}. The subject identifier SUBJECT and trait BMI. Since SUBJECT is not in the first column, `--phenotype_id SUBJECT` is required. Note that subject PGM in this example [phenotype file](files/phenotype.txt){:target="_blank"} is missing the BMI measurement.

<p><a id="covariates" title="Covariates" class="toc-item"></a></p>

### Covariates & Covariate Files

Covariate files are optional comma delimited files with header line with covariate name. Covariate files also assume the first column contains the individual ID unless the `--phenotype_id` is used. MMAP first searches the phenotype file for specified covariates, then the covariates files in order given. Thus, the first instance of the covariate is used if present in multiple files.

`--covariates <cov1> <cov2> ... <covN>` A list of covariates to be included in the analysis. If a covariate is not found in the phenotype file and/or covariate files, then the program exits. Individuals missing any covariate value are excluded from the analysis.

`--covariate_filename <file1> < file2> … <fileN>` Covariates specified by the `--covariates` option will be searched for in first the phenotype_file then sequentially through the list of covariate files. Thus, this option is only needed if covariates are not in the phenotype_file. The first file containing the covariate is used, thus there is no merging of within-covariate information across files. Covariate files can contain covariates with string values (so do not need to remove), but only covariates coded as numerical are currently supported. Missing data is coded as blank, so appear as “,,” in the csv file. Future options will include expanding a covariate into categorical values; for example, a single column containing four seasons would be expanded as 3 covariates. If no covariate files are given, then MMAP searches in the file specified by `--phenotype_filename` . If a covariate is not found, MMAP exits.

Example [covariate file](files/covariate.txt){:target="_blank"} contains data for 4 subjects from the pedigree file. The individual identifier SUBJECT must be the same as in the phenotype file. Subject 11 is missing AGE covariate.

NOTES:
- Caution when using Excel: If creating csv files form Excel it is recommended that a dummy column is added at the end because if there an individual has consecutive missing values in the file up the last column, sometimes the following rows get truncated. Excel will convert blanks to zeros when converting values. For example, if converting height from cm to inch by dividing the cm column by 2.54, missing cm values get a 0 in the inch column. Make sure no column has multiple values that are space/tab delimited as the column entries in the row will be shifted when exporting as csv.
- When moving a file from a Windows PC, run `dos2unix <windows-file>` to remove carriage returns.
- When moving a file from a MAC, run `mac2unix <mac-file>` to convert carriage returns to line feeds.

<p><a id="genotypes" title="Genotype Files" class="toc-item"></a></p>

### Genotype Files

MMAP uses binary only files for genotype analysis. For genome-wide association analysis the marker-by-subject (MxS) format is most efficient. MMAP provides utilities to convert comma separated text files with header to binary. The basic MxS format is as follows:
```
  SNPNAME            marker identifier
  RSNUM              rsnumber (or second SNP identifier)
  CHR	             1-22, X, Y, XY, MT
  POS	             position in base pairs
  STRAND	     +,- or blank
  NON-CODED_ALLELE   homozygote has dosage 0
  EFFECT_ALLELE	     homozygote has dosage 2
```
Only SNPNAME is required to be present, but default values of 0 for CHR and POS, ? for STRAND and 1 and 2 for the NON-CODED_ALLELE and EFFECT_ALLELE, respectively, will be entered if the tokens are not found. The file can contain any number of additional columns but must come before the list of individual IDs. These additional columns can be included in the genotype file and be used as filters or annotation by referencing the token in the header line in the appropriate command. Genotypes are stored in a variety of formats. Observed data is stored using 16 codes that represent phased and unphased states, partially typed states, and missing values. Thus, observed data is stored in a single byte. Imputed dosages can be stored as one or two bytes by scaling the value or as a double. One byte stores multiply the dosage by 100 so get 2 decimal place accuracy; two bytes multiplies by 10,000 so 4 place accuracy. If the original data is imputed dosage and not observed genotypes then the appropriate command line must be added to tell MMAP what data to expect. MMAP assumes comma-separated files (csv), but accepts space delimited also by adding `--genotype_space_delimiter` to the command line. MMAP accepts gzipped genotype files as input also. No additional command lines are necessary.

Example [genotype file](files/genotypes.txt){:target="_blank"} contains data for 5 subjects from the pedigree file using Affymetrix SNP ids as SNPNAME. The token STRAND is missing so the value will be set to ?. Marker SNP_A-2236359 has no rsnumber and is a deletion. Missing genotype values are coded as 3. There are two additional columns GENE and GROUP that can be included in the binary genotype file if desired. The following two commands are required to convert the text file to binary.

`--write_binary_genotype_file --csv_input_filename <input file> --binary_output_filename <output file>` converts the input file into MMAP binary format in the output file.

`--num_skip_fields` The number of columns to skip before the first subject ID. This command is required.

`--genotype_dosage_short` stores dosage as 10000*value, so suitable for dosage with 4 decimal place accuracy

`--genotype_dosage_char` stores dosage as 100*value, so suitable for dosage with 2 decimal place accuracy

`--genotype_dosage_double` stores dosage as value in file

`--genotype_space_delimiter`  add if genotype file is space rather than comma delimited

`--additional_marker_attributes <token1> <type1> … <tokenN> <typeN>` type is C for character string, D is double, I is integer. This option will include the additional columns from the genotype text file, if present.

`--output_marker_attribute <token1> <token2> … <tokenN>` The output file for GWAS will contain the standard tokens. Additional columns can be included using this option.

<p><a id="marker_options" title="Marker Options" class="toc-item"></a></p>

#### Options for marker analysis

`--chromosome <chrA> ... <chrN>` Analysis is restricted to chromosomes listed in the command line. Currently limited to autosomes Non-autosome chromosomes are designated by standard nomenclature: X, Y, XY, and MT.

`--genomic_region <chrA> <bp startA> <bp_endB> … <chrN> <bp startN> <bp_endN>` Analysis is restricted to genomic regions specified by chromosome and bp window. Base pair values. –chromosome 4 would be the same as –genomic_region 4 0 5000000.

`--marker_set <text file>` Analysis is restricted to the set of markers in the marker file. MMAP searches the SNPNAME and RSNUM columns to match the markers listed. NO header in the file


<p><a id="analysis" title="Analysis Set" class="toc-item"></a></p>

### Analysis set

MMAP takes the intersection of subjects in the pedigree and phenotype and also covariate and genotype files, if present, to generate the set of subjects used for analysis. Subjects with missing phenotype or covariate values are dropped. MMAP has an option to specify a subject set file, that if present, will also be included in the intersection. In the phenotype and covariate examples above the analysis set would for BMI with covariates AGE would be subject F from pedigree 1 and subject 12 from pedigree DM. If the genotype file is included then the analysis set contains only subject F.

`--subject_set <input file>` Single column file with no header that will control the individuals included in the analysis. This set is intersected with individuals with data from the phenotype and covariate file


<p><a id="running" title="Running MMAP" class="toc-item"></a></p>

### Running MMAP

MMAP requires that the binary relationship matrix be computed before any phenotype or genotype analysis. This matrix is then read in for other analyses.

`--compute_binary_relationship_matrix_by_groups <output file> --group_size <value>` Computes the pedigree specific relationship matrix (twice kinship) and stores in binary format to be read during analysis. Pre-computing this matrix is required to avoid recomputation for each analysis. The algorithm computes the matrix by groups to handle memory requirements of large pedigrees, as the matrix requires NxNx8 bytes of memory, where N is the number of individual. The binary file size is order NxNx4 bytes, thus depending on the application, it may be more efficient to restrict to the calculation to phenotyped individuals rather than the full pedigree by adding the subject_set option below. For example, when analyzing a Holstein pedigree of size ~240K with 60K genotyped animals only the 60K x 60K matrix was generated. For human pedigrees it is recommended to generate the full matrix as they rarely reach this size, even when the single_pedigree option is used.

Thus run the following command. Default group size is 1000.   
`mmap --ped <pedfile> --compute_binary_relationship_matrix_by_groups --binary_output-filname <output file> --group_size <value>`


<p><a id="output" title="Output Options" class="toc-item"></a></p>

### Output Options and Files:

`--file_suffix <string>` Adds string to output files to prevent clobbering from different analyses in the same directory

`--all_output` Generates additional output files that contain the likelihood values over the of h2 values, transformed phenotypes files.

`<trait>.<file_suffix>.poly.cov.csv` Contains trait statistics, number of observations, h2 estimate and p-value, beta, standard error and p-value of the fixed effects (uses t-test), percent variation the fixed effects account for, estimates of the total variance, additive variance, and error variance with standard errors. Also included is the standard errors of h2 estimate and standard error of the standard_deviation estimate (square root of variance estimate). MMAP uses the expected values of the information matrix, so no covariances between fixed and random effects are generated. Computing the matrix at the MLE estimates may be added in the future.

`<trait>.<file_suffix>.poly.cor.csv` The correlation/covariance between the beta estimators using the Fisher information matrix. The diagonal if the matrix is the standard deviation. Since expected values are used in the information matrix the covariance between fixed and random effects are assumed zero.

`<trait>.<file_suffix>.poly.model.csv` Contains the individuals used in the analysis, the observed phenotype, covariate values, fitted value and error residual calculated at the maximum likelihood estimate of h2. Other columns can be ignored and may be deleted in future versions. The column ERROR_RESIDUAL represents the residuals adjusted for the polygenic effect. These residuals can be treated as independent for analysis in programs that handle population samples.

`<trait>.<file_suffix>.spor.cov.csv` Same as poly version, but h2=0, so pedigree structure is ignored.

`<trait>.<file_suffix>.spor.cov.csv` Same as the poly version but h2=0, so pedigree structure is ignored.

<p><a id="example" title="Examples" class="toc-item"></a></p>

### Example Commands

`mmap --ped <pedfile> --read_binary_covariance_file pedigree.bin --trait HDL --phenotype_id MYEGO --phenotype_file pheno.csv --covariates AGE SEX BMI --covariate_file covarA.csv covarB.csv --file_suffix BMI`  
This command will analyze the trait HDL adjusting for covariates AGE, SEX and BMI. MMAP will look for the trait in pheno.csv then the covariates in the files pheno.csv, covarA.csv covarB.csv. The subject ID is assumed to by MYEGO in all three files. The output files with start with HDL and include BMI in the filename. The relationship matrix will be read from the file pedigree.bin

`mmap --ped <pedfile --read_binary_covariance_file pedigree.bin --trait HDL --phenotype_id MYEGO --phenotype_filename pheno.csv --covariates AGE SEX --file_suffix NO.BMI --binary_genotype_filename gwas.bin --model add --chromosome 4 X`  
This example is similar to that above except MMAP expects the covariates to be present in pheno.csv and BMI is dropped as a covariate. The file suffix is NO.BMI which will prevent the previous analysis results from being clobbered. MMAP will also perform marker analysis for chromosomes 4 and X using the additive model.

<p><a id="advanced" title="Advanced Topics" class="toc-item"></a></p>

| Advanced&nbsp;Topics&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | Documentation file |
| --- | --- |
| Working with the binary genotype file | [MMAP.genotype.pdf](files/MMAP.genotype.pdf){:target="_blank"} |
| Importing data from other formats: PLINK, MINIMAC, IMPUTE2, VCF | [MMAP.import.export.pdf](files/MMAP.import.export.pdf){:target="_blank"} |
| Genomic Relationship Matrices and PCs | [MMAP.genomic.matrix.pdf](files/MMAP.genomic.matrix.pdf){:target="_blank"} |
| Variance Component Estimation | [MMAP.variance.components.pdf](files/MMAP.variance.components.pdf){:target="_blank"} |
| Pedigree and Environmental Covariance Matrices | [MMAP.covariance.matrix.pdf](files/MMAP.covariance.matrix.pdf){:target="_blank"} |
| GxG, GxE, ExE Interaction Analysis and Sandwich Estimators | [MMAP.interaction.pdf](files/MMAP.interaction.pdf){:target="_blank"} |

<p><a id="genotypes_pdf" title="1. Genotype Files" class="toc-item"></a></p>

### MMAP Genotype File (from PDF)

MMAP has commands to manipulate the binary genotype file.
To transpose the file from marker-by-subject (MxS) to subject-by-marker (SxM) or subject-bymarker
(SxM) to marker-by-subject (MxS)  
`mmap --transpose_binary_genotype_file --binary_input_filename <infile> --binary_output_filename <outfile>`

If you transpose twice then you will get the original file.  
&nbsp; \- MxS format is useful for GWA, computation of LD matrices and allele frequency calculations.  
&nbsp; \- SxM format is useful for computation of genetic covariance matrices and haplotype analysis.

Creating a marker and/or subject reduced binary genotype file. Same command applies to MxS
and SxM formats. The outfile is the new binary genotype file.  
`mmap --write_reduced_genotype_binary --binary_input_filename <infile> --binary_output_filename <outfile> [marker and subject options]`

Example options:
```
--autosome               Extract the autosomal SNPs  
--chromosome <numbers>   Extract SNP on chromosomes in <numbers>  
--genomic_region <chr> <start bp> <stop bp>    Extract SNPs in the genomic region(s) specified.
--marker_set <file>      Extract markers in <file>
--subject_set <file>     Extract subjects in <file>
```

Extracting Genotypes as csv file. The infile must be MxS format. The outfile has the results. The options above are valid with this
command.  
`mmap --marker_by_subject_binary2csv --binary_input_filename <infile> --csv_output_filename <outfile>`

Allele frequency calculations.  The allele frequency for each SNP will be in the csv file.  
`mmap --write_binary_allele_frequency_file --binary_input_filename <MxS file> --binary_output_filename <filename> --csv_output_filename <filename>`


<p><a id="import" title="2. Genotype Import" class="toc-item"></a></p>

### MMAP Import Options (from PDF)

MMAP has commands to import data from Plink, Minimac, IMPUTE2 directly into a binary
genotype file and commands to export to Mach and Beagle format. Some of these options are
being beta tested or under development.

#### Plink

MMAP imports Plink binary format files into an SxM or MxS genotype binary file, depending on
the Plink format, which is automatically detected. Example below converts files \<prefix\>.bim,
\<prefix\>.bed, \<prefix\>.fam into binary genotype file \<mmap prefix\>.bin and MMAP
pedigree \<mmap prefix\>.ped.csv extracted from the \<prefix\>.fam.  
`mmap --plink_bfile2mmap –plink_bfile <prefix> -- binary_output_prefix <mmap prefix>`

#### Mach/MiniMac

MMAP imports Mach info and dosage files into an SxM binary genotype file. Since the map
information is not contained in the info file, the Mach map file is required. Dosage files can be
compressed or uncompressed. Options for reading in the probability file to create dominant and
recessive dosages are under development.  
`mmap --mach_dose2mmap –mach_info_filename <info file> --mach_dose_filename <dose file>
imputation_map_filename <map file> --binary_output_filename <SxM binary gentype
file> --genotype_dosage_short`  
\<info file> and \<dose file> are a list of Mach output files to be used. The chromosomes must be
in the same order in both files. The option `--genotype_dosage_short` stores the dosage as 2
bytes with precision 4-5 decimal places. The `--genotype_dosage_char` option will store the
dosage as 1 byte with 2-3 decimal place precision, reducing file size by half. It is recommended
to create a single binary file containing all the chromosomes for flexibility of analysis even
though the file will be large.

#### IMPUTE2

The IMPUTE2 import assumes that the probability .imputed and the information .imputed_info
files are available. Since subject id information is not contained in the output files, this
information must be included in the command line. This option is now set up to combine files per
chromosome to manage large files. Thus, the chromosome is required input. The default coding
is to use 1 byte. Since the probabilities are 3 decimal places, the 2 byte option is recommended.
The following conventions are used to handle the different variant types found in the current
imputation panels.  

Conventions:
1. SNPNAME is coded using the RSNUM value. If the marker was typed (2 in the output
column), the SNPNAME is the same as RSNUM, otherwise an “i” is appended. So rs123
becomes irs123
2. RSNUM is coded using the rs_id value in the .imputed or .imputed_info file, expect if it is
missing (dot in the output), then it is coded as <chr>:<position>. Use RSNUM when
reporting results.
3. STRAND is set to + by default as no strand information is available.
4. ALLELE as a coded as single characters using the nucleotides is both alleles are single
characters. Otherwise, R is used for the non-coded allele and I or D as the effect allele
depending on if the non-coded allele is a substring of the effect allele (I), or a superstring
(D). Marker with \<DEL\> are coded as R/R. MMAP outputs a file that contains the original
alleles and the codding. MMAP also supports multiallelic options where the alleles are as
the original but truncated beyond a maximum length.  
The following are required options:  
```
--impute2_prob2mmap
--impute2_prob_filename <file1> <file2> … <fileN>
--impute2_info_filename <file1> <file2> …<fileN>
--chromosome <chr>
--subject_id_filename <file>
--binary_output_filename <file>
--csv_output_filename <file>
```

The prob and info files should be in the same chunk order. The prob files can be gzip’d. The
subject id file does NOT have a header. The imputation quality score info is embedded in the
binary genotype file and can be outputted when running the single variant analysis. The default
is to create an additive dosage.

`--genotype_dosage_short` add this option to increase accuracy of stored dosage. Doubles file
size.  

Additional options to create alternative
dosage files:
```
--dominant_dosage creates dosage= 0*prob(AA) + 1*prob(AB) + 1*prob(BB)
--recessive_dosage creates dosage= 1*prob(AA) + 1*prob(AB) + 0*prob(BB)
--het_dosage creates dosage= 0*prob(AA) + 1*prob(AB) + 0*prob(BB)
```
**NOTE:** If you run `-–dom` with the binary_genotype_file created using `–dominant_dosage` you will
get a message that the option is not supported. For imputed data there is actually no model
since all genotypes have a value. The `-–dom` for observed data tells MMAP how to combine
genotypes into dominant dosages and to fill in missing values. These do not apply for in this
case. Thus, no model statement is needed.

Once the chromosome-specific files are created, they can then be combined into a single
MMAP binary using:
```
--combine_binary_genotype_files <file1> <file2> … <fileN>
--binary_output_filename <file>
```
The input files and output file are MMAP marker-by-subject binary genotype files.

**Recommendation:** To reduce the size of the combined binary genotype file, once the
chromosome specific files are created, run the allele frequency option. The output will contain
the minor allele frequency and imputation quality score, which can be used to extract a marker
set based on minor allele frequency and/or imputation quality threshold. This marker set can 
then be used to create reduced binary genotype files before combining to the full file. See
MMAP.genotype.pdf for details on these options.

#### VCF

Under development

<p><a id="export" title="3. Genotype Export" class="toc-item"></a></p>

### MMAP Export Options (from PDF)

#### Plink

`--subject_by_marker_mmap2plink --binary_input_filename <SxM binary gentype file> --
plink_output_prefix <prefix>`  
Creates \<prefix\>.map and \<prefix\>.ped which can then be converted into Plink binary format.
Currently no support of export directly into binary.

`--marker_by_subject_mmap2tped --binary_input_filename <SxM binary genotype file> --
plink_output_prefix <prefix>`  
Creates \<prefix\>.fam, \<prefix\>.bim and \<prefix\>.tped

#### Mach/Merlin

To be documented

#### Beagle

To be documented

#### MSMS

To be documented

#### ForSim

To be documented

#### Idcoeffs

To be documented

<br>

<p><a id="score_tests" title="Score Tests" class="toc-item"></a></p>

| Score&nbsp;Tests&nbsp;(including&nbsp;SKAT)&nbsp;&nbsp;&nbsp; | Documentation file |
|---|---|
| MMAP documentation for score tests | [MMAP.score.tests.pdf](files/MMAP.score.tests.pdf){:target="_blank"} |
| R script to convert MMAP prepScores output into a skatCohort object | [mmap2seqMeta.R](files/mmap2seqMeta.R) |
| MMAP scripts and examples for running score tests | [MMAP.score.test.tar.gz](files/MMAP.score.test.tar.gz) |
| MMAP snpinfo file used in CHARGE exome chip analysis | [<small>SNPInfo_HumanExome_12v1_rev5_AnalysisCols_noDups.tab</small>](files/SNPInfo_HumanExome_12v1_rev5_AnalysisCols_noDups.tab) |

<br>

<p><a id="cheatsheet" title="MMAP Cheat Sheet" class="toc-item"></a></p>

### MMAP Cheat Sheet

 &nbsp; &nbsp; &nbsp; &nbsp; [Click here for the MMAP Cheat Sheet](/CHEATSHEET.md){:target="_blank"}   
 &nbsp;
        
---

<p align="center">MMAP: Mixed Model Analysis for Pedigrees and Populations - Copyright © 2017</p>
&nbsp;

<!-- And now for the javascript... -->
  <script type="text/javascript" src="https://ajax.googleapis.com/ajax/libs/jquery/1.7.1/jquery.min.js"></script>
  <script type="text/javascript" src="/assets/js_custom/script.js"></script>
  <script type="text/javascript" src="/assets/js_custom/application.js"></script>
  
