% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bfoot.R
\name{combineTwoHotspots}
\alias{combineTwoHotspots}
\title{Combine two hotspot (peak) files by}
\usage{
combineTwoHotspots(csvfile1, csvfile2, name1, name2)
}
\arguments{
\item{bamfile}{Path to the input bam file}

\item{refgenome}{Reference genome, eg) "mm9", "hg19"}
}
\value{
Bedgraph filename
}
\description{
This function combines two hotspot files in the CSV format into
taking unions of hotspots
}
