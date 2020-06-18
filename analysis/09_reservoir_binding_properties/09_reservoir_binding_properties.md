Loading in peak\_occuence\_df and peak\_occurence\_matrix
=========================================================

``` r
base_path <- "../01_consensus_peaks/results"
peak_occurrence_matrix <- read.table(file.path(base_path,
                        "lncrna_mrna_promoter_peak_occurence_matrix.tsv"))

# Read in reservoir annotations combined with promoter peak occurence df. 
base_path <- "../08_defining_reservoirs/results"
peak_occurrence_df <- read.csv(file.path(base_path,
                                "08_peak_occurrence_df_promoter_types.csv"))
```

We are setting up to look at the binding distribution on expressed lncRNA and mRNA vs reservoir mRNA lncRNA binding. Distribution of reservoir vs all subsetted by lncRNA and mRNA

``` r
#### FIGURE: Figure 4C
g <- ggplot(peak_occurrence_df, aes(x = number_of_dbp))
g + geom_histogram(binwidth = 5)  + 
  xlim(30, 100) +
  facet_wrap(expression~gene_type, scales = "free_y")
```

    ## Warning: Removed 26802 rows containing non-finite values (stat_bin).

    ## Warning: Removed 8 rows containing missing values (geom_bar).

![](09_reservoir_binding_properties_files/figure-markdown_github/unnamed-chunk-2-1.png)

``` r
ggsave("figures/many_binders_histogram.png")
```

    ## Saving 7 x 5 in image

    ## Warning: Removed 26802 rows containing non-finite values (stat_bin).

    ## Warning: Removed 8 rows containing missing values (geom_bar).

``` r
ggsave("figures/many_binders_histogram.pdf")
```

    ## Saving 7 x 5 in image

    ## Warning: Removed 26802 rows containing non-finite values (stat_bin).

    ## Warning: Removed 8 rows containing missing values (geom_bar).

``` r
g <- ggplot(peak_occurrence_df %>% filter(number_of_dbp > 7), 
            aes(x = number_of_dbp, 
                fill = factor(reservoir, labels = c("Non-res", "Reservoir"))))
g + geom_density(alpha = 0.2) + 
  scale_fill_manual(values = c("#424242","#a8404c"),
                    name = " ") + 
  facet_grid(~gene_type) +
  ggtitle("DBP binding density",
          subtitle = "Promoters w/ > 7 DBPs")
```

![](09_reservoir_binding_properties_files/figure-markdown_github/unnamed-chunk-2-2.png)

``` r
ggsave("figures/many_binders_density_plot.png")
```

    ## Saving 7 x 5 in image

``` r
ggsave("figures/many_binders_density_plot.pdf")
```

    ## Saving 7 x 5 in image

We find that expressed or not expressed promoters still have similar distributions of DBP binding events. With some non-expressed promoters containing 100 independent DBP localization events.

Chi-Squared test for DBPs enriched and depleted at reservoirs versus non-reservoir promoters.

``` r
# Subset columns to only the reservoirs
res_gene_ids <- peak_occurrence_df$gene_id[peak_occurrence_df$reservoir == T]
res_occurrence_matrix <- peak_occurrence_matrix[ , res_gene_ids]
# Non reservoirs that also have more than 7 peaks bound
non_res <- peak_occurrence_df$gene_id[peak_occurrence_df$reservoir == F & 
                                        peak_occurrence_df$number_of_dbp > 7]
non_res_occurrence_matrix <- peak_occurrence_matrix[ , non_res]


# Number of promoters bound in each category
dbp_binding <- data.frame(dbp = rownames(peak_occurrence_matrix),
                          "num_res_promoters_bound" = rowSums(res_occurrence_matrix),
                          "num_nonres_promoters_bound" = rowSums(non_res_occurrence_matrix))

dbp_binding <- dbp_binding %>%
  mutate("num_res_promoters_not_bound" = ncol(res_occurrence_matrix) - num_res_promoters_bound,
         "num_nonres_promoters_not_bound" = ncol(non_res_occurrence_matrix) - num_nonres_promoters_bound)



for (i in 1:nrow(dbp_binding)) {
  df1 <- data.frame("gene_type" = c("reservoir","reservoir", "non_res", "non_res"),
                    "promoter_bound" = c("bound", "not_bound", "bound", "not_bound"),
                    "count" = c(dbp_binding$num_res_promoters_bound[i],
                                dbp_binding$num_res_promoters_not_bound[i],
                                dbp_binding$num_nonres_promoters_bound[i],
                                dbp_binding$num_nonres_promoters_not_bound[i])) %>%
    pivot_wider(names_from = gene_type, values_from = count) %>%
    column_to_rownames("promoter_bound") %>%
    as.matrix()
  
  csres <- chisq.test(df1)
  phi_coef <- phi(df1)
  
  # Add the results to the data frame
  dbp_binding[i, "chisq_stat"] <- csres$statistic
  dbp_binding[i, "chisq_pval"] <- csres$p.value
  dbp_binding[i, "reservoir_peaks_expected"] <- csres$expected["bound", "reservoir"]
  dbp_binding[i, "phi_coefficient"] <- phi_coef$phi
}
```

    ## Warning in chisq.test(df1): Chi-squared approximation may be incorrect

    ## Warning in stats::chisq.test(x, y, ...): Chi-squared approximation may be incorrect

    ## Warning in chisq.test(df1): Chi-squared approximation may be incorrect

    ## Warning in stats::chisq.test(x, y, ...): Chi-squared approximation may be incorrect

``` r
# Adjusting P value with BH correction, observed - expected and 
# writing this to a .csv file
dbp_binding$padj <- p.adjust(dbp_binding$chisq_pval, method = "BH")
dbp_binding$diff <- log2(dbp_binding$num_res_promoters_bound/dbp_binding$reservoir_peaks_expected)
write_csv(dbp_binding, "results/reservoir_chi_squared_results.csv")


# plotting Chi-squared test results 
g <- ggplot(dbp_binding, aes(x = diff, y = -log10(padj), label = dbp))
g + geom_point() +
  geom_hline(yintercept = -log10(0.001), lty = 2) +
  geom_vline(xintercept = -1, lty = 2) +
  geom_vline(xintercept = 1, lty = 2) +
  geom_vline(xintercept = 0, lty = 1) +
  geom_text_repel(data = dbp_binding %>% filter(-log10(padj) > 10 &
                    diff < -1)) + 
  ggtitle("Reservoir DBP bias",
          subtitle = "Chi-squared test res vs. non-res")
```

![](09_reservoir_binding_properties_files/figure-markdown_github/unnamed-chunk-3-1.png)

``` r
ggsave("figures/reservoir_dbp_chisq.pdf")
```

    ## Saving 7 x 5 in image

``` r
ggsave("figures/reservoir_dbp_chisq.png")
```

    ## Saving 7 x 5 in image

``` r
g <- ggplot(dbp_binding, aes(x = diff, y = -log10(padj), label = dbp))
g + geom_point() +
  geom_hline(yintercept = -log10(0.001), lty = 2) +
  geom_vline(xintercept = -1, lty = 2) +
  geom_vline(xintercept = 1, lty = 2) +
  geom_vline(xintercept = 0, lty = 1) +
  geom_text_repel(data = dbp_binding %>% filter(-log10(padj) > 3 & diff > 0)) + 
  ylim(0,5)
```

    ## Warning: Removed 89 rows containing missing values (geom_point).

    ## Warning: Removed 2 rows containing missing values (geom_text_repel).

![](09_reservoir_binding_properties_files/figure-markdown_github/unnamed-chunk-3-2.png)

``` r
ggsave("figures/reservoir_dbp_chisq_enriched.pdf")
```

    ## Saving 7 x 5 in image

    ## Warning: Removed 89 rows containing missing values (geom_point).

    ## Warning: Removed 2 rows containing missing values (geom_text_repel).

``` r
ggsave("figures/reservoir_dbp_chisq_enriched.png")
```

    ## Saving 7 x 5 in image

    ## Warning: Removed 89 rows containing missing values (geom_point).

    ## Warning: Removed 2 rows containing missing values (geom_text_repel).

``` r
# TODO ? UMAP reservoir promoter annotations ?
```
