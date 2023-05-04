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
  #     f(eGRN) = x(0.5y)
  
}
