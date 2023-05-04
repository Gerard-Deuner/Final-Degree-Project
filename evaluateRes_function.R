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

  # or 
  
  #         gene/TF         peak          resolution      validated            setting           dataset         corr.method     
  # link1   FOX1    chr1:187183-187634        0.5           1           combined_spearman     combined          spearman  
  # At the moment it the function is implemented for this format
  
  
  # Scoring function: (0-1)
  #     f(eGRN) = x(0.5y + 0.1z + 0.4t)
  
  
  # get number of validated TF-peak links from ChiP-seq
  ChiP.val <- val.df %>% 
    dplyr::filter(validation == "ChiP-seq") %>% # want to measure frequency of recoveries
    dplyr::select(resolution, setting, validated) %>%  # remove useless columns
    group_by(setting, resolution) %>% # group the data by resolution and setting
    summarise(recoveredTFPEAKs = sum(validated)) %>% # count number of TFs recovered per each resolution and setting
    dplyr::filter(resolution == res, setting == "combined_spearman_nomicro") %>% # filter for resolution and setting
    pull(recoveredTFPEAKs) # get number of TF-peak links

  # get number of validated peak-gene links from pcHi-C
  pcHiC.val <- val.df %>% 
    dplyr::filter(validation == "pcHi-C") %>% # want to measure frequency of recoveries
    dplyr::select(resolution, setting, validated) %>%  # remove useless columns
    group_by(setting, resolution) %>% # group the data by resolution and setting
    summarise(recoveredLINKs = sum(validated)) %>% # count number of TFs recovered per each resolution and setting
    dplyr::filter(resolution == res, setting == "combined_spearman_nomicro") %>% # filter for resolution and setting
    pull(recoveredLINKs) # get number of TF-peak links
  
  # get number of validated peak-gene links from eQTL
  eQTL.val <- val.df %>% 
    dplyr::filter(validation == "eQTL") %>% # want to measure frequency of recoveries
    dplyr::select(resolution, setting, validated) %>%  # remove useless columns
    group_by(setting, resolution) %>% # group the data by resolution and setting
    summarise(recoveredLINKs = sum(validated)) %>% # count number of TFs recovered per each resolution and setting
    dplyr::filter(resolution == res, setting == "combined_spearman_nomicro") %>% # filter for resolution and setting
    pull(recoveredLINKs) # get number of TF-peak links
  
  # get total number of inferred TF-peak links
  tf.peak.links <- val.df %>%                 
    dplyr::filter(validation == "ChiP-seq") %>% # want to measure frequency of recoveries 
    group_by(setting, resolution) %>% # group the data by resolution and setting
    summarise(recovered = sum(validated), GRN_links = length(validated)) %>% # count number of peak-gene links recovered per each resolution and setting
    dplyr::filter(resolution == res, setting == "combined_spearman_nomicro") %>% # filter for resolution and setting
    pull(GRN_links) # get number of total TF-peak links
  
  # get total number of inferred peak-gene links
  peak.gene.links <- val.df %>%                 
    dplyr::filter(validation == "pcHi-C") %>% # want to measure frequency of recoveries 
    group_by(setting, resolution) %>% # group the data by resolution and setting
    summarise(recovered = sum(validated), GRN_links = length(validated)) %>% # count number of peak-gene links recovered per each resolution and setting
    dplyr::filter(resolution == res, setting == "combined_spearman_nomicro") %>% # filter for resolution and setting
    pull(GRN_links) # get number of total peak-gene links
  
  # peak-gene QC score
  x <- qc
  
  # pcHi-C validation score
  y <- pcHiC.val / peak.gene.links
  
  # ChiP-seq validation score
  z <- ChiP.val / tf.peak.links
  
  # eQTL validation score
  t <- eQTL.val / peak.gene.links
  
  # compute function
  out <- x * ( (pcHiC.weight * y) + (ChiP.weight * z) + (eQTL.weight * t) )
  
  print(paste("Res", res, "Evaluation Score:", out, sep = " "))
}
