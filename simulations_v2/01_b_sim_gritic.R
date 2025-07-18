#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
n_clocks <- as.double(args[1])
n_events <- as.double(args[2])
pi = purity <- as.double(args[3])
coverage <- as.double(args[4])
epsilon <- as.double(args[5])
tolerance <- as.double(args[6])
max_attempts <- as.double(args[7])
seed <- as.double(args[8])
n_mutations <- as.numeric(args[9])
set.seed(seed)

source("utils.R")
source("utils_GRITIC.R")

MIN_MUTATIONS = 10 # write that is the minimum requested by gritic (check if it is like this)

print (paste0("n_clocks: ",n_clocks)) 
print (paste0("n_events: ",n_events))
print(paste0("purity: ",purity)) 
print(paste0("coverage: ",coverage))
print(paste0("epsilon: ",epsilon))
print(paste0("tolerance: ",tolerance)) 
print(paste0("max_attempts: ",max_attempts)) 
print(paste0("seed: ",seed))


library("reticulate")
library("dplyr")
use_python("/u/cdslab/scocomello/miniconda3/envs/r-reticulate/bin/python", required = TRUE)
py_config()
example_gritic <- reticulate::import_from_path("run_gritic_module", path = "../../")

safe_run <- function(expr, name) {
  tryCatch(
    expr,
    error = function(e) {
      msg <- paste(Sys.time(), "-", name, "failed with error:", e$message, "\n")
      writeLines(msg, error_log)
      return(NULL)
    }
  )
}

main_dir = paste0("results_gritic/sim_", n_clocks, "_", n_events, "_", purity, "_", coverage, "_", n_mutations, "/")
if (!dir.exists(main_dir)) {
  dir.create(main_dir)  
}

for (i.iter in 1:5) {
  sub_dir = paste0(main_dir, i.iter)

  error_log <- file(paste0(sub_dir, "/error_log_2.txt"), open = "wt")
  
  sim = readRDS(paste0(sub_dir,"/sim.rds"))
  
  N_events = nrow(sim$cn)
  cn <- sim$cn %>% 
    dplyr::rename(Major=nMaj1_A, minor=nMin1_A, from=startpos, to=endpos)
  muts <- sim$muts %>% 
    dplyr::rename(from=start, to=end)
  input_data_CNAqc <- CNAqc::init(mutations=muts, 
                                  cna = cn, 
                                  purity = pi, 
                                  ref = "hg19")
  x = list(
    cna = cn, 
    mutations = muts,
    metadata = dplyr::tibble(sample = "sample", purity=pi)
  )
  # input_data = list (mutations=muts, cna = cn, purity=pi)
  data_GRITIC = convert_CNAqc_to_GRITIC(input_data_CNAqc) 
  mutation_table = paste0(sub_dir, "/mutation_tsv.tsv")
  copy_number_table = paste0(sub_dir, "/cna_tsv.tsv")
  write_tsv(data_GRITIC$mut, mutation_table)
  write_tsv(data_GRITIC$cn, copy_number_table)
  
  output_dir = paste0(sub_dir, "/results_gritic/")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)  
  }
  safe_run(fit_gritic(mutation_path=mutation_table, copy_number_path = copy_number_table,purity = pi, sample_id = "S123",
                      output_dir = paste0(".", output_dir),
                      wgd_status = TRUE,
                      penalty = FALSE,
                      subclone_path = NULL,
                      sex = "XY",
                      merge_cn=FALSE,
                      plot_trees = TRUE), "fit_gritic")
  
  res_gritic <- read.csv(paste0(sub_dir,"/results_gritic/S123_gain_timing_table_wgd_segments.tsv"), sep = "\t")
  
  # Save simulation results
  if (!is.null(res_gritic)) saveRDS(res_gritic, paste0(sub_dir, "/res_gritic.rds"))
  
  if (!is.null(res_gritic)) merged_res = readRDS(paste0(sub_dir, "/merged_res.rds"))
  # add gritic results
  if (!is.null(res_gritic)) merged_res = merged_res %>% mutate("tau_gritic"=res_gritic$Timing) #check that the correspondance of the segments are correct, by hand they are
  if (!is.null(res_gritic)) saveRDS(merged_res, paste0(sub_dir, "/merged_res.rds"))
  
  # Close the error log file
  close(error_log)

}
