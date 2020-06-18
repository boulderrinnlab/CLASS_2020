Reservoir nascent transcription
===============================

Goal: To analyze nascent RNA-sequencing data by PRO-SEQ methods at reservoir promoters. We are reasoning that althought there is not a "mature transcript" produced from these promoters that "nascent transcription = transcripts being produced by PolR2 but may not finalize to a "mature transcript"

First we will count the "nascent" reads in the promoter regions defined across this study (note: nascent transcripts are often found in promoters and in equal amounts on + and - strand of DNA)

``` r
#Defining promoter transcriptional start sites from gencode_gr
gencode_gr <- rtracklayer::import("/Shares/rinn_class/data/genomes/human/gencode/v32/gencode.v32.annotation.gtf")

# Organizing to get' 'all' promoters (lncRNA and mRNA) +/- 3Kb 
# from start site (from get_promoter_regions function)
promoters_df <- rtracklayer::import("../01_consensus_peaks/results/lncrna_mrna_promoters.gtf") %>%
  as.data.frame() %>%
  dplyr::select("gene_id", "seqnames", "start", "end", "strand") 
names(promoters_df) <- c("GeneID", "Chr", "Start", "End", "Strand")
promoters_df[which(promoters_df$Start < 0), "Start"] <- 0
write.table(promoters_df, "results/gencode_v32_promoters_6kb.SAF", sep = "\t",
            quote = F, row.names = F)
```

Using Rsubread to quantify the nascent epxression over promoter windows (3kb +/- start) for analysis in 11\_reservoir\_chromatin\_properties.

``` r
# Rsubread to make counts on features.
counts <- featureCounts(c("data/516.sorted.bam", "data/829.sorted.bam"), 
                          annot.ext = "results/gencode_v32_promoters_6kb.SAF", 
                          isPairedEnd = TRUE, nthreads = 16)
```

    ## 
    ##         ==========     _____ _    _ ____  _____  ______          _____  
    ##         =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
    ##           =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
    ##             ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
    ##               ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
    ##         ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
    ##        Rsubread 2.0.1
    ## 
    ## //========================== featureCounts setting ===========================\\
    ## ||                                                                            ||
    ## ||             Input files : 2 BAM files                                      ||
    ## ||                           o 516.sorted.bam                                 ||
    ## ||                           o 829.sorted.bam                                 ||
    ## ||                                                                            ||
    ## ||              Annotation : gencode_v32_promoters_6kb.SAF (SAF)              ||
    ## ||      Dir for temp files : .                                                ||
    ## ||                 Threads : 16                                               ||
    ## ||                   Level : meta-feature level                               ||
    ## ||              Paired-end : yes                                              ||
    ## ||      Multimapping reads : counted                                          ||
    ## || Multi-overlapping reads : not counted                                      ||
    ## ||   Min overlapping bases : 1                                                ||
    ## ||                                                                            ||
    ## ||          Chimeric reads : counted                                          ||
    ## ||        Both ends mapped : not required                                     ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## //================================= Running ==================================\\
    ## ||                                                                            ||
    ## || Load annotation file gencode_v32_promoters_6kb.SAF ...                     ||
    ## ||    Features : 36814                                                        ||
    ## ||    Meta-features : 36814                                                   ||
    ## ||    Chromosomes/contigs : 25                                                ||
    ## ||                                                                            ||
    ## || Process BAM file 516.sorted.bam...                                         ||
    ## ||    WARNING: Single-end reads were found.                                   ||
    ## ||    Total alignments : 14965166                                             ||
    ## ||    Successfully assigned alignments : 2136291 (14.3%)                      ||
    ## ||    Running time : 0.04 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file 829.sorted.bam...                                         ||
    ## ||    WARNING: Single-end reads were found.                                   ||
    ## ||    Total alignments : 27527321                                             ||
    ## ||    Successfully assigned alignments : 4154920 (15.1%)                      ||
    ## ||    Running time : 0.05 minutes                                             ||
    ## ||                                                                            ||
    ## || Write the final count table.                                               ||
    ## || Write the read assignment summary.                                         ||
    ## ||                                                                            ||
    ## \\============================================================================//

``` r
# Converting to Reads Per Kilobase of gene length.
rpk <- counts$counts / (counts$annotation$Length/1000)

# Making an expression summary data frame 
# rpk_m > first RPK is converted to rpk_M (rpk per millions of reads run)
# TPM > second rpk_m is converted to TPM (Transcripts per millions of reads run )

expression <- data.frame("rpk" = rpk) %>%
  rownames_to_column("gene_id") %>%
  pivot_longer(2:3, names_to = "sample", values_to = "rpk")
expression_summary <- expression %>% 
  group_by(sample) %>%
  summarize(total_rpk = sum(rpk, na.rm = T))
```

    ## `summarise()` ungrouping output (override with `.groups` argument)

``` r
expression_summary$rpk_m <- expression_summary$total_rpk / 1e6
expression <- merge(expression, expression_summary)
expression$tpm <- expression$rpk / expression$rpk_m
tpm <- expression %>% group_by(gene_id) %>%
  summarize(tpm = mean(tpm, na.rm = T))
```

    ## `summarise()` ungrouping output (override with `.groups` argument)

``` r
genes <- gencode_gr[gencode_gr$type == "gene"] %>% as.data.frame() %>%
  dplyr::select(gene_name, gene_id)
tpm <- merge(tpm, genes)
tpm <- tpm %>% dplyr::select(gene_id, gene_name, tpm)
write_csv(tpm, "results/k562_nascent_tpm.csv")
```

``` r
# Let's now add nascent_tpm to our 'peak_occurrence_df' saving as "peak_occurrence_df_nascent"
# First loading the 'peak_occurence_df' data frame we created in 08_defining reservoirs

peak_occurrence_df_08 <- read.csv("../08_defining_reservoirs/results/08_peak_occurrence_df_promoter_types.csv")
peak_occurrence_nascent <- read.csv("results/k562_nascent_tpm.csv")

# Adding in nascent by gene_id (and gene name) -- double checked k562_nascent_tpm vs output and is all good --
peak_occurrence_df_08$nascentTpm <- peak_occurrence_nascent$tpm

# Writting out .csv of peak_occurrence
write_csv(peak_occurrence_df_08, "results/peak_occurrence_df_nascent_added.csv")

# This .csv will be read into 11_reservoir_chromatin_properties.
```
