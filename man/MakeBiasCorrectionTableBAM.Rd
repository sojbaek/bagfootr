% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bfoot.R
\name{MakeBiasCorrectionTableBAM}
\alias{MakeBiasCorrectionTableBAM}
\title{Generate n-mer bias tables from a given bam file}
\usage{
MakeBiasCorrectionTableBAM(
  bamfile = "",
  outfile = "",
  refgenome = "",
  np = 6,
  mapdir = "",
  atac = F
)
}
\arguments{
\item{bamfile}{Path to the input bam file}

\item{outfile}{Output Filename}

\item{refgenome}{Reference genome, eg) "mm9", "hg19"}

\item{np}{Number of base pairs, 2 = dimer, 4 = tetramer, 6 = hexamer (default)}

\item{mapdir}{Directory to the mappability files}

\item{atac}{Logical, T for ATAC-seq data}
}
\value{
Output Filename
}
\description{
This function generate n-mer bias tables from a given BAM file.
If the mappability files directory is given, mappability is used to
generate the bias table.  The Reference genome must be given to run this function.
}
