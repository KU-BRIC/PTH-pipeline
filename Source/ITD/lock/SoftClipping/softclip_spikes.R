suppressMessages({
	library(tidyverse)
    library(zoo)
})


options(tibble.width = Inf)


###
### Command Line Arguments
###

args = commandArgs(trailingOnly = TRUE)

sample_name <- args[1L]
out_dir <- args[2L]
sample_chr_pos_ref_cigar_seq_path <- args[3L]   # tsv with sample pos cigar sequence
MIN_CIGAR_SPIKE <- as.integer(args[4L])           # min cigar spike size for it to be considered a potential ITD
REF_GEN <- (args[5L])
read_length <- as.integer(args[6L])


cat(
  sprintf("
          arg1 = %s\n
          arg2 = %s\n
          arg3 = %s\n
          arg4 = %s\n
          arg5 = %s\n
	  arg6 = %s\n",
	  args[1L], args[2L], args[3L], args[4L], args[5L], args[6L])
)




###
### Function Definitions
###

my_rollsuml <- function(v, window) {
  #v <- 1:20
  #window <- 4
  sum_v <- vector(mode = "integer", length = length(v))

  for (i in 1:20) {
    sum_v[i] <- sum(v[i:(i+window-1)], na.rm = TRUE)
  }

  return(sum_v)
}



process_sample_pos_cig_seq <- function(dat_cpcsf, read_length) { #sample, chr, pos, ref, cigar, seq, .pb = NULL) {
  # TBD: add avg non-0 cigar height
  #
 print("process_sample_pos_cig_seq")
 print(read_length) 
 print(class(read_length))
 wind <- as.numeric(read_length) * 2

  dat_cpcsf %>% #head(100000) %>%
    filter(cigar != "*") %>%
    
    # if it starts with eg 36S and the string contains no Ns, get the number into S_n
    mutate(S_n = map2_int(cigar, seq,
                          function(cg, s)
                            ifelse(str_detect(cg, "^[:digit:]+S") & str_detect(s, "N+", negate = TRUE),
                                   as.integer(str_extract(cg, "[:digit:]+")),
                                   as.integer(0)))) %>%
    
    # if S_n > 0, get the substring
    mutate(S_seq = map2_chr(S_n, seq,
                            function(n, sq)
                              ifelse(n > 0, str_sub(sq, 1, n), ""))) %>%

    group_by(pos, sample) %>%
    mutate(ctypes = str_c(unique(cigar), collapse = ":"),
           n_ctypes = length(cigar),
           n_unique_ctypes = length(unique(cigar)),
           n_reads = length(pos),
           n_softclipped = sum(map_lgl(cigar, ~str_detect(., "^[:digit:]+S"))),
           n_match = sum(map_lgl(cigar, ~str_detect(., "^[:digit:]+M$")))) %>%

    ungroup %>%
    group_by(sample, chr, pos, n_match, n_reads) %>%
    nest %>%
    ungroup %>%
    mutate(n_2RL = rollapply(n_reads, wind, sum, align="center", fill = 0)) %>%
    unnest(data) %>%
    mutate(vaf_estim = round(n_softclipped/n_2RL, 3))
}


mark_cigar_spikes <- function(dat, MCS) {
  # adds boolean feature 'max_cigar' that 
  cat("min_cigar_spike = ", MCS, "\n")
  cat(".......", class(MCS), " ...... ")
  dat %>% head %>% print
  dat %>%
    group_by(pos, sample) %>%
    mutate(max_cigar = n_unique_ctypes >= MCS & n_unique_ctypes == max(n_unique_ctypes)) %>%
    ungroup %>%
    unique
}



plot_sample <- function(dat, fid, seq_u = "", seq_d = "", seq_ins = "", VAF = "", title = "") {
  max_y <- dat %>% filter(sample == fid, max_cigar) %>% dplyr::select(n_unique_ctypes) %>% unique %>% unlist %>% max
  pp <- dat %>%
    filter(sample == fid) %>%
    ggplot(aes(x = pos, y = n_unique_ctypes)) +
    geom_col(width = 3, position = position_dodge())+
    geom_text(aes(x=pos, y=n_unique_ctypes+2, label = ifelse(max_cigar, pos, "")), size = 3) +
    geom_text(aes(x=pos, y=n_unique_ctypes+1, label = ifelse(max_cigar, n_unique_ctypes, "")), size = 3) +
    #geom_text(aes(x=pos, y=n_unique_ctypes+6, label = ifelse(max_cigar, str_c("NORMAL: ", seq_u, "-", seq_d, sep = ""), "")), size = 3) +
    #geom_text(aes(x=pos, y=n_unique_ctypes+5, label = ifelse(max_cigar, str_c("SPLICE: ", seq_ins, "-", seq_d, sep = ""), "")), color = "red", size = 3) +
    #geom_text(aes(x=pos, y=n_unique_ctypes+4, label = ifelse(max_cigar, str_c("VAF: ", VAF), "")), color = "red", size = 3) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, size = 5),
      plot.margin = unit(rep(0.5, 4), "cm")
    ) +
    ylim(0, max_y + 6) +
    ggtitle(title)
  
  # print needed else it doesn't hv time to plot to pdf
  return(pp)
}



###
### Look for SoftClipped Spikes
###

sample_chr_pos_ref_cigar_seq <- read_tsv(sample_chr_pos_ref_cigar_seq_path, col_names = c("sample", "chr", "pos", "ref", "alt", "cigar", "seq"), col_types = "ccicccc") %>%
  replace_na(replace = list(seq = ""))
  sample_chr_pos_ref_cigar_seq %>% head %>% print
dat_tmp <- process_sample_pos_cig_seq(sample_chr_pos_ref_cigar_seq, read_length)
dat_tmp %>% head %>% print
cat("BBB\n")
dat_processed <- mark_cigar_spikes(dat = dat_tmp, MCS = MIN_CIGAR_SPIKE)
cat("CCC\n")
cigar_spike_samples <- dat_processed %>% #filter(sample == "Horizon-CMP016") %>%
  dplyr::filter(max_cigar) %>% dplyr::select(sample, chr, pos, ref, alt, n_2RL, n_softclipped, vaf_estim, n_unique_ctypes, max_cigar) %>% unique
cat("DDD\n")
print(nrow(cigar_spike_samples))
cigar_spike_samples <- cigar_spike_samples# %>% mutate(pos = pos - 1)

###
### Plot and Write Results
###


if (nrow(cigar_spike_samples) > 0) {
  out_pdf <- file.path(out_dir, str_c(sample_name, ".pdf"))
  cat(out_pdf)
  pdf(out_pdf)
  print(plot_sample(dat_processed, dat_processed$sample, title = dat_processed$sample))
  dev.off()
 
  out_tsv <- file.path(out_dir, str_c(sample_name, ".tsv"))
  cat(out_tsv)
  write_tsv(cigar_spike_samples, file = out_tsv)
}


