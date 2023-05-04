################################################
# Score Resolutions based on previous analyses #
################################################

# create function to evaluate which is the best resolution based on validation scores
evaluateRes <- function(data,                  # df with resolutions as rows and validation stats as columns. (data.frame)
                        res,                   # resolution to evaluate. (int)
                        qc,                    # parameter that indicates if that resolution has passed the peak-gene QC filters. 0 = No. 1 = Yes. (int)
                        pcHiC.weight = 0.5,    # weight given to the pcHi-C validation. (int). Default = 0.5. 
                        ChiP.weight = 0.1,     # weight given to the ChiP-seq validation. (int). Default = 0.1.
                        eQTL.weight = 0.4      # weight given to the eQTL validation. (int). Default = 0.4 
) {
  
  # data format:
  #   res  pcHiC.val   ChiP.val    eQTL.val    peak.gene.links    tf.peak.links
  #   0.1     130         3           85          1320                  754
  #   ...     ...        ...         ...          ...                   ...

  # Scoring function: (0-1)
  #     f(eGRN) = x(0.5y + 0.1z + 0.4t)
  
  # peak-gene QC score
  x = qc
  
  # pcHi-C validation score
  y = data[res, pcHiC.val] / data[res, peak.gene.links]
  
  # ChiP-seq validation score
  z = data[res, ChiP.val] / data[res, tf.peak.links]
  
  # eQTL validation score
  t = data[res, eQTL.val] / data[res, peak.gene.links]
  
  # compute function
  out <- x * ( (pcHiC.weight * y) + (ChiP.weight * z) + (eQTL.weight * t) )
  
  print(paste(res, "Evaluation Score:", out, sep = " "))
}
