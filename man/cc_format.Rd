% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/models.R
\name{CC_format}
\alias{CC_format}
\title{Format data for use by CopyClust function}
\usage{
CC_format(data_input, reference_genome = "hg18", probes = 1e+05)
}
\arguments{
\item{data_input}{A data frame generated from the output of \href{https://bioconductor.org/packages/release/bioc/html/DNAcopy.html}{DNAcopy}. The data frame must contain six columns with the following column names: \code{ID}, \code{chrom}, \code{loc.start}, \code{loc.end}, \code{num.mark}, and \code{seg.mean}.}

\item{reference_genome}{Parameter for reference genome. Formats the genomic ranges used as features for \code{\link[=CopyClust]{CopyClust()}} to the appropriate reference genome. Valid inputs are \code{hg18}, \code{hg19}, and \code{hg38}. Default is \code{hg18}.}

\item{probes}{Parameter for number of probes to use to calculate feature values. Specified number of probes will be used to calculate values for model features selected from all available probes equally spaced across the genome. Default is \code{100000} probes. A greater number of probes decreases the processing speed. Providing a values greater than the number of available probes will results in a error.}
}
\value{
A data frame with sample IDs as rows and 478 model features as columns that can be used with the \code{\link[=CopyClust]{CopyClust()}} function.
}
\description{
Formats raw data from DNACopy format into 478 genomic range features required to run the \code{\link[=CopyClust]{CopyClust()}}.
Reference genome (\code{hg18}, \code{h19}, or \code{hg38}) must be specified with the \code{reference_genome} parameter.
}
\examples{
data("test_data_raw")
data_for_CopyClust = CC_format(test_data_raw, reference_genome = "hg19", probes = 100000)

}
\seealso{
\code{\link[=CopyClust]{CopyClust()}}, \href{https://bioconductor.org/packages/release/bioc/html/DNAcopy.html}{DNAcopy}
}
\author{
Cameron C. Young
}
