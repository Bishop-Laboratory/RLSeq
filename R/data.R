#' Auxiliary Data
#'
#' A list containing data used by RLSeq functions. It can also be useful
#' for checking the available modes and genomes in RLSeq. See also the
#' `data-raw/auxdata.R` script that was used to create it.
#'
#' ## Structure
#'
#' A named list containing the following entries:
#'
#' * `db_cols`
#'   - A `tbl` with colors associated with each database in RLHub
#'   useful for plotting. See also [RLHub::annotations].
#' * `annotypes`
#'   - A `tbl` containing the annotation databases and annotation types
#'   available from RLBase. See also [RLHub::annotations].
#' * `ip_cols`
#'   - A `tbl` containing the colors associated with each "Immunoprecipitation
#'   type" (ip_type) in RLBase. See also [RLHub::rlbase_samples].
#' * `mode_cols`
#'   - A `tbl` containing the colors associated with each R-loop mapping
#'   mode in [RLHub::rlbase_samples].
#' * `heat_cols`
#'   - A `tbl` containing the colors associated with user-supplied data and
#'   RLBase data when running [corrHeatmap].
#' * `label_cols`
#'   - A `tbl` containing the colors associated with the labels in
#'   [RLBase](https://gccri.bishop-lab.uthscsa.edu/rlbase/).
#'   See also [RLHub::rlbase_samples].
#' * `prediction_cols`
#'   - A `tbl` containing the colors associated with the predictions in
#'   [RLBase](https://gccri.bishop-lab.uthscsa.edu/rlbase/).
#'   See also [RLHub::rlbase_samples].
#' * `prediction_label_cols`
#'   - A `tbl` containing the colors associated with the prediction-label
#'   combinations in
#'   [RLBase](https://gccri.bishop-lab.uthscsa.edu/rlbase/).
#'   See also [RLHub::rlbase_samples].
#' * `available_modes`
#'   - A `tbl` containing the modes available in RLBase and associated metadata.
#'   See also [RLHub::rlbase_samples].
#' * `available_genomes`
#'   - A `character` showing all the official UCSC genomes available for use
#'   with RLSeq. See also [available_genomes].
#' * `misc_modes`
#'   - A `character` showing the R-loop mapping modes that are lumped into
#'   the 'misc' category for simplification of plotting.
#'
#' @examples
#' auxdata
#'
#' @export
"auxdata"

#' Available Genomes
#'
#' Contains metadata about all the genomes available in UCSC. It contains
#' derived metadata, such as the effective genome sizes as well. See also
#' the `data-raw/available_genomes.R` script to see processing steps.
#'
#' ## Structure
#'
#' `available_genomes` is a `data.frame` with the following columns:
#'
#' * `UCSC_orgID`
#'   - Official UCSC ID of the genome
#' * `description`
#'   - Verbose description of the assembly, source, and year/month of entry.
#' * `nibPath`
#'   - Endpoint of the genome in UCSC gbdb.
#' * `organism`
#'   - Name of the organism.
#' * `defaultPos`
#'   - Default location of genome browser view for this genome.
#' * `active`
#'   - Description not available.
#' * `orderKey`
#'   - Description not available.
#' * `genome`
#'   - The name of the genome.
#' * `scientificName`
#'   - The scientific name of the organism.
#' * `htmlPath`
#'   - Path in UCSC gbdb to the `description.html` file for the genome.
#' * `hgNearOk`
#'   - Description not available.
#' * `hgPbOk`
#'   - Description not available.
#' * `sourceName`
#'   - Name of organization providing the genome.
#' * `taxId`
#'   - The taxonomy ID of the organism.
#' * `genes_available`
#'   - If `TRUE`, the gene annotations are available in GTF format.
#' * `year`
#'   - The year the genome assembly was added.
#' * `eff_genome_size_XXbp`
#'   - The effective genome size of this genome. Calculated at various read
#'   lengths with [khmer](https://khmer.readthedocs.io/en/latest/)
#'   and used to improve the accuracy of analysis. See
#'   the `data-raw/available_genomes.R` script to see how this calculation
#'   was performed.
#' * `genome_length`
#'   - The total length of the genome.
#' * `rlfs_available`
#'   - If `TRUE`, R-loop forming sequences annotations are available in the
#'   RLBase AWS S3 repository.
#'
#' @examples
#' available_genomes
#'
#' @export
"available_genomes"


#' Genome Masks
#'
#' A collection of genome masks for use with [analyzeRLFS]. See the
#' `data-raw/genome_masks.R` script for the processing steps.
#'
#' ## Structure
#'
#' `genomeMasks` is a named list of `GRanges` objects. Each entry in the
#' list follows the naming convention: `<genome>.masked`, where `<genome>`
#' is an official UCSC genome ID. Each entry contains
#' a `GRanges` object with the masked ranges from `<genome>`. The genomes
#' provided correspond to the masked genomes available in
#' [BSgenome::available.genomes].
#'
#' @examples
#' genomeMasks
#'
#' @export
"genomeMasks"


#' Random Genomic Windows
#'
#' A collection of random genomic windows for use with [noiseAnalyze]. See the
#' `data-raw/genome_masks.R` script for the processing steps.
#'
#' ## Structure
#'
#' `randomWindows` is a named list of `tbl` objects containing ~1000 random
#' genomic regions. The names are the genomes to which the regions correspond.
#'
#' Columns (mirrors
#' [bed3 format](https://genome.ucsc.edu/FAQ/FAQformat.html#format1)).
#'
#' 1. `chrom` - the chromosome in UCSC format
#' 2. `start` - the starting position
#' 3. `end` - the end position
#'
#' @examples
#' randomWindows
#'
#' @export
"randomWindows"


#' RLBase Sample Transcript Feature Overlaps
#'
#' Summary statistics from transcript feature overlap analysis of peaks
#' from all RLBase samples.
#'
#' ## Structure
#'
#' `rlsampleTxOl` is a `tbl` with the following columns:
#'
#' - `rlsample`
#'   + The RLBase sample identifier for the sample.
#'   + Matches the `rlsample` column in [RLHub::rlbase_samples]
#' - `feature`
#'   + The transcript feature for which overlap analysis was performed.
#'   + These features were derived from the Transcript Features collection
#'   described in [RLHub::annotations]
#' - `n`
#'   + The raw number of peaks from the sample overlapping a feature.
#' - `pct`
#'   + The proportion of peaks fro the sample overlapping a feature.
#'
#' @examples
#' rlsampleTxOl
#'
#' @export
"rlsampleTxOl"


#' RLBase Sample Noise Analysis Results
#'
#' Average signal from noise analysis of RLBase samples. See [noiseAnalyze] for
#' more detail regarding how signal was initially calculated.
#'
#' ## Structure
#'
#' `rlbaseNoiseAnalyze` is a `tbl` with the following columns:
#'
#' - `rlsample`
#'   + The RLBase sample identifier for the sample.
#'   + Matches the `rlsample` column in [RLHub::rlbase_samples]
#' - `value`
#'   + The mean signal from noise analysis. See [noiseAnalyze].
#'
#' @examples
#' rlbaseNoiseAnalyze
#'
#' @export
"rlbaseNoiseAnalyze"
