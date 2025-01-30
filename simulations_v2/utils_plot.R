
ALPHA = .8
color = c(
  'AmplificationTimeR' = alpha('forestgreen', alpha = ALPHA),
  'MutationTimeR' = alpha('steelblue', alpha = ALPHA),
  'tickTack' = alpha('orange', alpha = ALPHA),
  'tickTackH' = alpha('firebrick3', alpha = ALPHA)
)


convert_name = function(n) {
  if (grepl("AmpTime", n)) return("AmplificationTimeR")
  if (grepl("MutTimeR", n)) return("MutationTimeR")
  if (grepl("tickTack_h", n)) return("tickTackH")
  if (grepl("tickTack", n)) return("tickTack")
  stop("error name not recognized")
}

# Percent error plots ####
plot_over_purity = function(all_res) {
  r = all_res %>%
    tidyr::pivot_longer(c(tau_AmpTimeR, tau_MutTimeR, tau_tickTack, tau_tickTack_h)) %>%
    dplyr::group_by(name, setting, purity) %>%
    #dplyr::summarise(metric = sqrt(mean(value - true_tau)**2)) %>% # RMSE
    dplyr::summarise(metric = mean(abs(value - true_tau) / true_tau), .groups = "drop") # PERCENT ERROR
  r$name = lapply(r$name, convert_name) %>% unlist()

  r %>%
    ggplot(mapping = aes(x = as.factor(purity), y=metric, fill=name)) +
    geom_boxplot() +
    theme_bw() +
    scale_y_continuous(transform = "log10") +
    labs(x = "Purity", y="Percent error", fill="") +
    scale_fill_manual(values = color)
}

plot_over_coverage = function(all_res) {
  r = all_res %>%
    tidyr::pivot_longer(c(tau_AmpTimeR, tau_MutTimeR, tau_tickTack, tau_tickTack_h)) %>%
    dplyr::group_by(name, setting, coverage) %>%
    #dplyr::summarise(metric = sqrt(mean(value - true_tau)**2)) %>% # RMSE
    dplyr::summarise(metric = mean(abs(value - true_tau) / true_tau), .groups = "drop") # PERCENT ERROR
  r$name = lapply(r$name, convert_name) %>% unlist()

  r %>%
    ggplot(mapping = aes(x = as.factor(coverage), y=metric, fill=name)) +
    geom_boxplot() +
    theme_bw() +
    scale_y_continuous(transform = "log10") +
    labs(x = "Coverage", y="Percent error", fill="") +
    scale_fill_manual(values = color)
}

plot_over_nmutations = function(all_res) {
  r = all_res %>%
    tidyr::pivot_longer(c(tau_AmpTimeR, tau_MutTimeR, tau_tickTack, tau_tickTack_h)) %>%
    dplyr::group_by(name, setting, n_mutations) %>%
    #dplyr::summarise(metric = sqrt(mean(value - true_tau)**2)) %>% # RMSE
    dplyr::summarise(metric = mean(abs(value - true_tau) / true_tau), .groups = "drop") # PERCENT ERROR
  r$name = lapply(r$name, convert_name) %>% unlist()

  r %>%
    ggplot(mapping = aes(x = as.factor(n_mutations), y=metric, fill=name)) +
    geom_boxplot() +
    theme_bw() +
    scale_y_continuous(transform = "log10") +
    labs(x = "N mutations", y="Percent error", fill="") +
    scale_fill_manual(values = color)
}

plot_nclocks_v_nevents = function(all_res) {
  r = all_res %>%
    tidyr::pivot_longer(c(tau_AmpTimeR, tau_MutTimeR, tau_tickTack, tau_tickTack_h)) %>%
    dplyr::group_by(name, setting, n_clocks, n_events) %>%
    #dplyr::summarise(metric = sqrt(mean(value - true_tau)**2)) %>% # RMSE
    dplyr::summarise(metric = mean(abs(value - true_tau) / true_tau), .groups = "drop") # PERCENT ERROR
  r$name = lapply(r$name, convert_name) %>% unlist()

  r %>%
    ggplot(mapping = aes(x = as.factor(n_events), y=metric, fill=name)) +
    geom_boxplot() +
    theme_bw() +
    ggh4x::facet_nested(~"N tau"+n_clocks) +
    scale_y_continuous(transform = "log10") +
    labs(x = "N segments", y="Percent error", fill="") +
    scale_fill_manual(values = color)
}

# Rand index plot ####
plot_rand_index = function(res_clust) {
  colnames(res_clust)[7:10] = c("AmplificationTimeR", "MutationTimeR", "tickTack", "tickTackH")
  res_clust %>%
    dplyr::select(n_clocks, n_events, AmplificationTimeR, MutationTimeR, tickTack, tickTackH) %>%
    tidyr::pivot_longer(c(AmplificationTimeR, MutationTimeR, tickTack, tickTackH)) %>%
    ggplot(mapping = aes(x = name, y=value, fill=name)) +
    geom_boxplot() +
    ggh4x::facet_nested("N segments"+n_events~"N tau"+n_clocks) +
    scale_y_continuous(breaks = scales::pretty_breaks(n=3)) +
    theme_bw() +
    scale_fill_manual(values = color) +
    labs(fill = "", x = "", y="Rand Index") +
    theme(
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank()
    )
}
