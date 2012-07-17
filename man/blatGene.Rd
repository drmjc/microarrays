\name{blatGene}
\alias{blatGene}
\title{Affy probe sequences to BLAT input}
\usage{
  blatGene(affyid, probe, filename)
}
\arguments{
  \item{affyid}{Affy probeset ID}

  \item{probe}{BioC probe package name}

  \item{filename}{output file name}
}
\value{
  a fasta file, suitable for pasting into the BLAT tool @
  UCSC UCSC BLAT:
  \url{http://genome.ucsc.edu/cgi-bin/hgBlat?command=start}
  Code written by Jim McDonald, on the bioconductor mailing
  list: Re: [BioC] Affy chip annotation changes, 3 November
  2010 12:08:14 AM AEDT
}
\description{
  Function to obtain the probe sequences for an Affy
  probeset ID & writes it out in a format suitable for
  UCSC's BLAT
}
\author{
  Mark Cowley
}

