#' Affy probe sequences to BLAT input
#' 
#' Function to obtain the probe sequences for an Affy probeset ID & writes it
#' out in a format suitable for UCSC's BLAT
#' 
#' @param affyid Affy probeset ID
#' @param probe BioC probe package name
#' @param filename output file name
#' @return a fasta file, suitable for pasting into the BLAT tool @@ UCSC UCSC
#'   BLAT: \url{http://genome.ucsc.edu/cgi-bin/hgBlat?command=start} Code written by
#'   Jim McDonald, on the bioconductor mailing list: Re: [BioC] Affy chip
#'   annotation changes, 3 November 2010 12:08:14 AM AEDT
#' @author Mark Cowley
#' @export
blatGene <- function(affyid, probe, filename){
   ## affyid == Affy probeset ID
   ## probe == BioC probe package name
   ## filename == output file name
   require(probe, quietly = TRUE, character.only = TRUE)
   tmp <- data.frame(get(probe))
   if(length(affyid) > 1){
       seqnc <- vector()
       for(i in seq(along = affyid))
           seqnc <- c(seqnc, tmp[tmp$Probe.Set.Name == affyid[i], 1])
   }else{
       seqnc <- tmp[tmp$Probe.Set.Name == affyid,1]
   }
   out <- vector()
   if(length(seqnc) > 25) warning("Blat will only return values for 25 or fewer sequences!",
                                  call. = FALSE)
   for(i in seq(along = seqnc)) out <- rbind(out, rbind(paste("> Probe", i, sep=""), seqnc[i]))
   write.table(out, filename, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}
