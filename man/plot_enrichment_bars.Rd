% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/function.R
\name{plot_enrichment_bars}
\alias{plot_enrichment_bars}
\title{Visualize enrichment results as a bar plot}
\usage{
plot_enrichment_bars(data, pvalue_threshold = 0.05)
}
\arguments{
\item{data}{A dataframe containing enrichment results with columns: ID, pvalue, and optionally cell_type.}

\item{module_info}{A dataframe mapping module names to cell types (required if cell_type not in data).
Should contain columns: name, cell_type.}
}
\value{
A ggplot object showing enrichment results visualization.
}
\description{
Creates a bar plot of enrichment results ordered by -log10(p-value) and colored by cell type.
The function automatically merges cell type information if not already present in the data.
}
