source("./simulate_functions.R")
starting='1:1'
karyotype_1='4:2' 
karyotype_2='4:2'
CCF_1=1
purity=.8
b_model = suppressWarnings(branching_evolution(starting, karyotype_1, karyotype_2, CCF_1, purity)) %>%
  Reduce(f = bind_rows)


first_child = karyotype_1
second_child = karyotype_2

l1_model = NULL
if ((grepl("0", karyotype_1) & grepl("0", karyotype_2)) | !grepl("0", karyotype_1)){
  l1_model = suppressWarnings(linear_evolution(starting, karyotype_1, karyotype_2, CCF_1, purity)) %>%
  Reduce(f = bind_rows)
}
