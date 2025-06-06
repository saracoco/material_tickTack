ALPHA = .8
# color = c(
#   'AmplificationTimeR' = alpha('forestgreen', alpha = ALPHA),
#   'MutationTimeR' = alpha('steelblue', alpha = ALPHA),
#   'tickTack' = alpha('orange', alpha = ALPHA),
#   'tickTackH' = alpha('firebrick3', alpha = ALPHA)
# )

color = c(
  'AmplificationTimeR' = alpha('#868686', alpha = ALPHA),
  'MutationTimeR' = alpha('#EFC000', alpha = ALPHA),
  'tickTack baseline' = alpha('#7AA6DC', alpha = ALPHA),
  'tickTack full' = alpha('#CD534C', alpha = ALPHA)
)

k_colors = list(
  '1:1' = '#228B22CC',
  '1:0' = 'steelblue',
  '2:0' = 'turquoise4',
  '2:1' = ggplot2::alpha('orange', .8),
  '2:2' = 'firebrick3'
)

convert_name = function(n) {
  if (grepl("AmpTime", n)) return("AmplificationTimeR")
  if (grepl("MutTimeR", n)) return("MutationTimeR")
  if (grepl("tickTack_h", n)) return("tickTack full")
  if (grepl("tickTack", n)) return("tickTack baseline")
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
    geom_boxplot(lwd=.3, outlier.size = .5) +
    theme_bw() +
    scale_y_continuous(transform = "log10") +
    labs(x = "Purity", y="Log percent error", fill="") +
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
    geom_boxplot(lwd=.3, outlier.size = .5) +
    theme_bw() +
    scale_y_continuous(transform = "log10") +
    labs(x = "Coverage", y="Log percent error", fill="") +
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
    dplyr::mutate(name = factor(name, levels=names(color))) %>%
    ggplot(mapping = aes(x = as.factor(n_mutations), y=metric, fill=name)) +
    geom_boxplot(lwd=.3, outlier.size = .5) +
    theme_bw() +
    scale_y_continuous(transform = "log10") +
    labs(x = "N mutations", y="Log percent error", fill="") +
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
    dplyr::mutate(name = factor(name, levels=names(color))) %>%
    ggplot(mapping = aes(x = as.factor(n_events), y=metric, fill=name)) +
    geom_boxplot(lwd=.3, outlier.size = .5) +
    theme_bw() +
    ggh4x::facet_nested(~"N tau"+n_clocks) +
    scale_y_continuous(transform = "log10") +
    labs(x = "N segments", y="Log percent error", fill="") +
    scale_fill_manual(values = color)
}

# Rand index plot ####
plot_rand_index = function(res_clust) {
  colnames(res_clust)[7:10] = c("AmplificationTimeR", "MutationTimeR", "tickTack baseline", "tickTack full")
  res_clust %>%
    dplyr::select(n_clocks, n_events, AmplificationTimeR, MutationTimeR, `tickTack baseline`, `tickTack full`) %>%
    tidyr::pivot_longer(c(AmplificationTimeR, MutationTimeR, `tickTack baseline`, `tickTack full`)) %>%
    dplyr::mutate(name = factor(name, levels=names(color))) %>%
    ggplot(mapping = aes(x = name, y=value, fill=name)) +
    geom_boxplot(lwd=.3, outlier.size = .5) +
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




get_reference = function(ref)
{
  if(ref %in% c("hg19", "GRCh37"))
    return(CNAqc::chr_coordinates_hg19)

  if(ref %in% c("hg38", "GRCh38"))
    return(CNAqc::chr_coordinates_GRCh38)

  stop("Available references: hg19 (or GRCh37), and hg38 (or GRCh38)")
}


add_breakpoints_to_plot = function(segments, base_plot, max_Y_height, circular)
{
  off_plot = segments %>% dplyr::filter(total > max_Y_height)
  if (nrow(off_plot) > 0)
  {
    base_plot = base_plot +
      ggplot2::geom_hline(
        yintercept = max_Y_height,
        size = .2,
        color = 'darkgray',
        linetype = 'dashed'
      )
    L = ggplot2::ggplot_build(base_plot)$layout$panel_params[[1]]
    Lx = abs(L$x.range[2] - L$x.range[1]) * .85

    base_plot = base_plot +
      ggplot2::geom_label(
        data = data.frame(
          x = Lx,
          y = L$y.range[2] - 0.5,
          label = paste0('< ', max_Y_height)
        ),
        ggplot2::aes(x = x, y = y, label = label),
        fill = 'darkgray',
        color = 'white',
        size = 2,
        nudge_y = 0,
        nudge_x = 0,
        inherit.aes = FALSE
      )
  }

  L = ggplot2::ggplot_build(base_plot)$layout$panel_params[[1]]

  if(!circular){
    if (L$y.range[2] <= 5) {
      base_plot = base_plot + ggplot2::ylim(-0.5, 5)
    }
  }

  # =-=-=-=-=-=-=-=-=-=-=-=-
  # Breakpoints annotations
  # =-=-=-=-=-=-=-=-=-=-=-=-
  breakpoints = data.frame(
    x = segments$from,
    y = 0.1,
    outern = segments$Major > max_Y_height
  )

  base_plot = base_plot +
    ggplot2::geom_point(
      data = breakpoints %>% filter(!outern),
      ggplot2::aes(x = x, y = y),
      size = .5,
      shape = 1,
      color = 'darkgray'
    ) +
    ggplot2::geom_point(
      data = breakpoints %>% filter(outern),
      aes(x = x, y = y),
      size = .5,
      color = 'black'
    )

  return(base_plot)
}

get_drivers = function(x,
                       chromosomes = paste0('chr', c(1:22, 'X', 'Y')),
                       which = 'VAF')
{
  if(!has_driver_data(x)) return(NULL)

  # CCF?
  if (all(is.null(x$CCF_estimates)) & which == "CCF")
  {
    warning("Input does not have CCF estimates, see ?compute_CCF to determine CCF values.")
    return(NULL)
  }

  drivers_list = NULL

  if(which == "VAF")
    drivers_list = x$mutations %>%
    dplyr::filter(is_driver, chr %in% chromosomes)

  if(which == "CCF")
    drivers_list = CNAqc::CCF(x) %>%
    dplyr::filter(is_driver, chr %in% chromosomes)

  return(drivers_list)
}


has_driver_data = function(x)
{

  cn = colnames(x$mutations)

  if (all(c("is_driver", "driver_label") %in% cn)) return(TRUE)

  if ("is_driver" %in% cn & !("driver_label" %in% cn)){
    cli::cli_warn("Column 'is_driver' is annotated but 'driver_label' no -- did you try to add drivers data?")
  }

  if ("driver_label" %in% cn & !("is_driver" %in% cn)){
    cli::cli_warn("Column 'driver_label' is annotated but 'is_driver' no -- did you try to add drivers data?")
  }

  return(FALSE)
}




add_drivers_to_segment_plot = function(x, drivers_list, base_plot)
{
  if(is.null(drivers_list)) return(base_plot)
  if (nrow(drivers_list) == 0) return(base_plot)

  L = ggplot2::ggplot_build(base_plot)$layout$panel_params[[1]]

  drivers_list = CNAqc:::relative_to_absolute_coordinates(
    x,
    drivers_list %>% dplyr::filter(is_driver)
  )

  drivers_list$y = L$y.range[2] * 0.9

  base_plot +
    ggplot2::geom_vline(
      data = drivers_list,
      show.legend = FALSE,
      ggplot2::aes(xintercept = from),
      linetype = 'dashed',
      color = 'black',
      size = .3
    ) +
    ggrepel::geom_label_repel(
      data = drivers_list,
      ggplot2::aes(
        x = from,
        y = y,
        label = driver_label,
      ),
      ylim = c(L$y.range[2] * .9, NA),
      size = 2,
      nudge_y = 0,
      nudge_x = 0,
      show.legend = FALSE
    ) +
    ggplot2::coord_cartesian(clip = 'off')
}




merge_timing_and_segments <- function(x, chromosomes = paste0('chr', c(1:22)), max_Y_height = 6, cn = 'absolute', highlight = x$most_prevalent_karyotype, highlight_QC = FALSE) {

  K = x$K

  cnaqc_x = CNAqc::init(mutations = x$mutations, cna = x$cna, purity = x$metadata$purity, ref = x$reference_genome)

  plot_CNA = plot_segments_tick_tack_CN(cnaqc_x, K = K) +
    theme(legend.position='right', panel.spacing = unit(0, "lines")) +
    labs(caption = NULL) +
    ggplot2::theme(axis.title.x = element_blank())

  x$reference_genome = "hg38"
  data_plot <- plot_segments_tick_tack_data(x, K = K) +
    theme(legend.position='right',panel.spacing = unit(0, "lines")) +
    ggplot2::theme(axis.title.x = element_blank())

  timing_plot <- plot_segments_tick_tack(x, colour_by = "clock_mean", K = K) +
    theme(legend.position='right',panel.spacing = unit(0, "lines"))

  vaf_plot <- plot_vaf(x, K) +
    theme(plot.background = element_rect(fill = "white", color = "white", size = 20))

  # segment_plot <- plot_segments_h(x, chromosomes, max_Y_height, cn, highlight, highlight_QC) +
  #   ggplot2::theme(axis.title.x = element_blank())  # Keep chromosome labels only on this plot

  pA = timing_plot + CNAqc:::my_ggplot_theme() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      legend.position = "left"
      )
  pB = plot_CNA + CNAqc:::my_ggplot_theme() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      legend.position = "left"
    )
  pC = data_plot + CNAqc:::my_ggplot_theme() +
    theme(legend.position = "left") +
    labs(y = "VAF")
  pD = vaf_plot + CNAqc:::my_ggplot_theme()

  des_left = "
  AAAA#
  AAAAE
  AAAAE
  AAAAE
  AAAAE
  BBBBE
  BBBBE
  BBBBE
  CCCCE
  CCCCE
  DDDDD"

  pp = pA + pB + pC + guide_area() + pD +
    plot_layout(design = des_left, guides = "collect") &
    theme(legend.position = "bottom", legend.direction = "horizontal")
  pp
}

plot_segments_tick_tack <- function(x, colour_by = "clock_mean", K = 1) {

  results <- x$results_timing
  reference_genome <- get_reference(x$reference_genome)
  segments <- results$data$accepted_cna
  segments$from = lapply(segments$segment_name, function(s) {unlist(strsplit(s, "_"))[2]}) %>% unlist() %>% as.numeric()
  segments$to = lapply(segments$segment_name, function(s) {unlist(strsplit(s, "_"))[3]}) %>% unlist() %>% as.numeric()

  vfrom = reference_genome$from
  names(vfrom) = reference_genome$chr

  absoulte_segments <- segments %>%
    dplyr::mutate(from = .data$from + vfrom[.data$chr],
                  to = .data$to + vfrom[.data$chr])

  summarized_results <- results$draws_and_summary[[K]]$summarized_results %>%
    dplyr::mutate(from = absoulte_segments[.data$segment_id,]$from) %>%
    dplyr::mutate(to = absoulte_segments[.data$segment_id,]$to) %>%
    dplyr::mutate(tau_mean = ifelse(.data$clock_mean  < 1, .data$clock_mean , 1)) %>%
    dplyr::mutate(tau_high = ifelse(.data$clock_high < 1, .data$clock_high, 1)) %>%
    dplyr::mutate(tau_low = .data$clock_low)
  ##############

  k_colors = list(
    '1:1' = '#228B22CC',
    '1:0' = 'steelblue',
    '2:0' = 'turquoise4',
    '2:1' = ggplot2::alpha('orange', .8),
    '2:2' = 'firebrick3'
  )

  # six_color_palette <- six_color_palette <- RColorBrewer::brewer.pal(6, "RdYlBu")
  summarized_results <- summarized_results %>%
    dplyr::mutate(tau_cluster = factor(
      .data[[colour_by]],
      levels = sort(unique(.data[[colour_by]])),
      labels = paste(seq_along(sort(unique(.data[[colour_by]]))))
    ))

  capt_label = paste0("Tumour type ", x$ttype,
                      " Ploidy ", x$ploidy, "; Purity  ", x$purity,
                      '; n = ', x$results_timing$data$input_data$N, ' accepted mutations in ',
                      x$results_timing$data$input_data$S,
                      ' accepted segments'
  )

  my_palette <- c(  "#66a61e",  "#7570b3", "#e7298a", "#1b9e77", "#d95f02")

  p <- CNAqc:::blank_genome('hg19',chromosomes = paste0("chr", c(1:22)))+
    ggplot2::geom_segment(
      data = summarized_results,
      ggplot2::aes(
        y = .data$clock_mean,
        yend = .data$clock_mean,
        x = .data$from,
        xend = .data$to
      )
    ) +
    ggplot2::geom_rect(
      data = summarized_results,
      ggplot2::aes(
        xmin = .data$from,
        xmax = .data$to,
        ymin = .data$tau_low,
        ymax = .data$tau_high,
        fill = .data$tau_cluster
      ),
      alpha = 0.5
    ) +
    ggplot2::scale_fill_manual(values = my_palette)+
    ggplot2::theme(
      legend.position = "bottom",
      axis.text.x = ggplot2::element_text(angle = 0)
    ) +
    ggplot2::ggtitle(label = "", subtitle = capt_label) +
    ggplot2::labs(
      x = "Chromosome",
      y = bquote("Pseudotime"~tau),
      fill = "Cluster"
    )

  drivers_list = get_drivers(x, chromosomes = paste0('chr', c(1:22)))
  base_plot = add_drivers_to_segment_plot(x, drivers_list = drivers_list, p)
  base_plot
}

plot_segments_tick_tack_CN <- function(cnaqc_x, K = K, max_alleles = 6, chromosomes = paste0("chr", c(1:22))) {

  reference_genome <- get_reference(cnaqc_x$reference_genome)

  vfrom = reference_genome$from
  names(vfrom) = reference_genome$chr

  segments <- cnaqc_x$cna %>%
    dplyr::filter(chr %in% chromosomes)

  absolute_segments <- segments %>%
    dplyr::mutate(from = .data$from + vfrom[.data$chr],
                  to = .data$to + vfrom[.data$chr])

  absolute_segments$karyotype = paste(absolute_segments$Major, absolute_segments$minor, sep = ":")
  absolute_segments = absolute_segments %>%
    dplyr::filter(Major <= max_alleles & minor <= max_alleles)

  k_colors = list(
    '1:1' = '#228B22CC',
    '1:0' = 'steelblue',
    '2:0' = 'turquoise4',
    '2:1' = ggplot2::alpha('orange', .8),
    '2:2' = 'firebrick3'
  )
  p = CNAqc:::blank_genome(cnaqc_x$reference_genome, chromosomes = paste0("chr", c(1:22)))

  # Add shadows
  p = p + ggplot2::geom_rect(
    data = absolute_segments %>% dplyr::filter(karyotype %in% c('2:0', '1:0', '1:1', '2:1', '2:2')),
    ggplot2::aes(
      xmin = from,
      xmax = to,
      ymin = -Inf,
      ymax = Inf,
      fill = factor(karyotype, levels = c('2:0', '1:0', '1:1', '2:1', '2:2'))
    ),
    alpha = .3
  ) +
    ggplot2::scale_fill_manual(values = k_colors) +
    ggplot2::guides(fill = ggplot2::guide_legend('', override.aes = list(alpha = 1)))

  # Add segments
  p = p +
    ggplot2::geom_segment(
      data = absolute_segments %>%
        dplyr::mutate(Major = as.numeric(Major) + .1, minor = as.numeric(minor) - .1) %>%
        dplyr::select(karyotype, chr, from, to, Major, minor) %>%
        tidyr::pivot_longer(!c(karyotype, chr, from, to)),
      ggplot2::aes(
        x = from,
        xend = to,
        y = value,
        color = as.factor(.data$name)
      ),
      size=1.5
    ) +
    ggplot2::scale_color_manual(values = c("Major" = "red", "minor" ="steelblue")) +
    ggplot2::guides(color = ggplot2::guide_legend('')) +
    labs(y = "Allel count") +
    ylim(c(0, max_alleles))

  p
}


plot_segments_tick_tack_data <- function(x, colour_by = "clock_mean", K = K) {

  mutations <- x$mutations
  results <- x$results_timing
  reference_genome <- get_reference(x$reference_genome)

  vfrom = reference_genome$from
  names(vfrom) = reference_genome$chr

  absolute_mutations <- mutations  %>%
    dplyr::mutate(from = .data$from + vfrom[.data$chr],
                  to = .data$to + vfrom[.data$chr])

  segments <- results$data$accepted_cna
  segments$from = lapply(segments$segment_name, function(s) {unlist(strsplit(s, "_"))[2]}) %>% unlist() %>% as.numeric()
  segments$to = lapply(segments$segment_name, function(s) {unlist(strsplit(s, "_"))[3]}) %>% unlist() %>% as.numeric()

  absolute_segments <- segments %>%
    dplyr::mutate(from = .data$from + vfrom[.data$chr],
                  to = .data$to + vfrom[.data$chr])

  summarized_results <- results$draws_and_summary[[K]]$summarized_results %>%
    dplyr::mutate(from = absolute_segments[.data$segment_id,]$from) %>%
    dplyr::mutate(to = absolute_segments[.data$segment_id,]$to) %>%
    dplyr::mutate(tau_mean = ifelse(.data$clock_mean  < 1, .data$clock_mean , 1)) %>%
    dplyr::mutate(tau_high = ifelse(.data$clock_high < 1, .data$clock_high, 1)) %>%
    dplyr::mutate(tau_low = .data$clock_low)

  accepted_mutations = data.frame()
  for (segment_idx in 1:nrow(summarized_results)) {
    segment <- summarized_results[segment_idx, ]
    print(segment$chr)
    segment_mutations <- absolute_mutations %>%
      dplyr::filter(.data$chr == segment$chr, .data$from > segment$from, .data$to < segment$to) %>%
      tidyr::drop_na(DP)
    print(nrow(segment_mutations))
    # if (nrow(segment_mutations)> 40){
    accepted_mutations <- bind_rows(accepted_mutations, segment_mutations)
    # }
  }

  matched_mutations <- accepted_mutations

  k_colors = list(
    '2:0' = 'turquoise4',
    '2:1' = ggplot2::alpha('orange', .8),
    '2:2' = 'firebrick3'
  )

  CNAqc:::blank_genome(ref=x$reference_genome, chromosomes = paste0("chr", c(1:22)))+
    ggplot2::geom_point(
      data = matched_mutations,
      ggplot2::aes(
        x = from,
        y = .data$NV / DP,
        color = as.factor(.data$karyotype)
      ),
      alpha = 0.5,
      size=0.1
    ) +
    ggplot2::scale_color_manual(values = k_colors, name = "CN") +
    ggplot2::guides(
      color = ggplot2::guide_legend(override.aes = list(size = 4))) +
    ggplot2::labs(
      x = "Chromosome",
      y = bquote("Variant Allele Frequency (VAF)")
    )
}



plot_segments_h = function(x,
                           chromosomes = paste0('chr', c(1:22)),
                           max_Y_height = 6,
                           cn = 'absolute',
                           highlight = x$most_prevalent_karyotype,
                           highlight_QC = FALSE,
                           ...)
{


  # Segments

  segments_original = x$cna %>%
    dplyr::filter(chr %in% chromosomes)  %>%
    dplyr::mutate(
      total = Major + minor,
      karyotype = paste0(Major, ':', minor)
    )

  segments_original = CNAqc:::relative_to_absolute_coordinates(x, segments_original)

  #
  # results <- x$results_timing
  # reference_genome <- get_reference (x$reference_genome)
  # vfrom = reference_genome$from
  # names(vfrom) = reference_genome$chr
  #
  # segments <- results$data$accepted_cna
  # segments$from = lapply(segments$segment_name, function(s) {unlist(strsplit(s, "_"))[2]}) %>% unlist() %>% as.numeric()
  # segments$to = lapply(segments$segment_name, function(s) {unlist(strsplit(s, "_"))[3]}) %>% unlist() %>% as.numeric()
  #
  # absolute_segments <- segments %>%
  #   dplyr::mutate(from = .data$from + vfrom[.data$chr],
  #                 to = .data$to + vfrom[.data$chr])
  #
  # summarized_results <- results$draws_and_summary[[K]]$summarized_results %>%
  #   dplyr::mutate(from = absolute_segments[.data$segment_id,]$from) %>%
  #   dplyr::mutate(to = absolute_segments[.data$segment_id,]$to) %>%
  #   dplyr::mutate(tau_mean = ifelse(.data$clock_mean  < 1, .data$clock_mean , 1)) %>%
  #   dplyr::mutate(tau_high = ifelse(.data$clock_high < 1, .data$clock_high, 1)) %>%
  #   dplyr::mutate(tau_low = .data$clock_low)
  #
  # segments_filter_total = data.frame()
  # for (segment_idx in 1:nrow(summarized_results)) {
  #   segment <- summarized_results[segment_idx, ]
  #   print(segment$chr)
  #   segments_filter <- segments_original %>%
  #     dplyr::filter(.data$chr == segment$chr, .data$from >= segment$from, .data$to <= segment$to) %>%
  #   print(nrow(segments_filter))
  #   # if (nrow(segment_mutations)> 40){
  #   segments_filter_total <- bind_rows(segments_filter_total, segments_filter)
  #   # }
  # }
  #

  # df <- segments_filter_total %>%
  #   mutate(chr_num = as.numeric(str_extract(chr, "\\d+")))

  # Find chr with max and min numeric values
  # chr_max <- df %>% filter(chr_num == max(chr_num)) %>% pull(chr_num) %>% unique()
  # chr_min <- df %>% filter(chr_num == min(chr_num)) %>% pull(chr_num) %>% unique()
  #
  # chromosomes = paste0('chr', c(chr_min:chr_max))

  # Standard plot -- baseline genome reference
  base_plot = blank_genome(chromosomes = chromosomes,
                           ref = x$reference_genome)
  segments_original = segments_original %>%
    dplyr::filter(chr %in% chromosomes)  %>%
    dplyr::mutate(
      total = Major + minor,
      karyotype = paste0(Major, ':', minor)
    )


  segments <- segments_original
  # =-=-=-=-=-=-=-=-=-=-=-=-
  # Shadow for highligthing
  # =-=-=-=-=-=-=-=-=-=-=-=-
  if(highlight_QC & !is.null(x$peaks_analysis)) {
    base_plot = CNAqc:::add_shadow_to_plot_QC(segments, base_plot)
  } else {
    base_plot = CNAqc:::add_shadow_to_plot(segments, base_plot, highlight)
  }



  # =-=-=-=-=-=-=-=-=-=-=-=-
  # Draw Segments
  # =-=-=-=-=-=-=-=-=-=-=-=-
  base_plot = add_segments_to_plot(
    segments = segments %>% dplyr::filter(total <= max_Y_height),
    base_plot = base_plot,
    cn = cn)

  # Extract subclonal segments
  subclonal_segments = NULL
  if (!is.null(x$cna_subclonal) & nrow(x$cna_subclonal) > 0)
  {
    subclonal_segments = x$cna_subclonal %>%
      dplyr::filter(chr %in% chromosomes)

    if (nrow(subclonal_segments) > 0)
    {
      base_plot = add_subclonal_segments_to_plot(
        segments = subclonal_segments %>%
          relative_to_absolute_coordinates(x = x),
        base_plot = base_plot,
        cn = cn
      )
    }
  }

  # Fragmentation ~ add some annotation to hihglight that
  if (!is.null(x$arm_fragmentation))
  {
    fragmented = x$arm_fragmentation$table %>%
      dplyr::filter(significant, chr %in% chromosomes) %>%
      dplyr::mutate(label = paste0(chr, arm)) %>%
      dplyr::pull(label)

    if (length(fragmented) > 0)
    {
      expanded_reference = CNAqc:::expand_reference_chr_to_arms(x) %>%
        dplyr::filter(chr %in% fragmented)

      # base_plot = base_plot +
      #   geom_rect(
      #     data = expanded_reference,
      #     aes(
      #       xmin = from,
      #       xmax = to,
      #       ymin = -Inf,
      #       ymax = Inf
      #     ),
      #     alpha = .2,
      #     fill = NA,
      #     color = 'purple4'
      #   )
      #
      # base_plot +
      #   geom_segment(
      #     data = expanded_reference,
      #     aes(
      #       x = from,
      #       xend = to,
      #       y = -0.2,
      #       yend = -0.2,
      #     ),
      #     color = 'purple4',
      #     linetype = 1,
      #     size = 2
      #   )

      base_plot = base_plot +
        ggplot2::geom_label(
          data = expanded_reference,
          ggplot2::aes(
            x = from,
            label = gsub(pattern = 'chr', replacement = '', chr),
            y = -0.2
          ),
          fill = 'purple4',
          color = 'white',
          linewidth = 2,
          label.padding = unit(0.05, 'cm')
        )
    }
  }

  # Add extras
  if (!is.null(x$arm_fragmentation))
    capt_label = paste0(
      capt_label,
      '; ',
      x$arm_fragmentation$table %>%
        dplyr::filter(significant, chr %in% chromosomes) %>%
        length,
      ' fragmented arms'
    )
  # =-=-=-=-=-=-=-=-=-=-=-=-
  # Breakpoints annotations
  # =-=-=-=-=-=-=-=-=-=-=-=-
  # base_plot = CNAqc:::add_breakpoints_to_plot(segments, base_plot, max_Y_height, circular = FALSE)

  # =-=-=-=-=-=-=-=-=-=-=-=-
  # Drivers annotations
  # =-=-=-=-=-=-=-=-=-=-=-=-
  drivers_list = get_drivers(fit, chromosomes = chromosomes)
  # if(!circular)
  base_plot = add_drivers_to_segment_plot(x,
                                          drivers_list = drivers_list,
                                          base_plot)



  return(base_plot)
}


blank_genome = function(ref = "GRCh38", chromosomes = paste0('chr', c(1:22, 'X', 'Y')), label_chr = -0.5, cex = 1){
  reference_coordinates = get_reference(ref) %>%
    filter(chr %in% chromosomes)

  low = min(reference_coordinates$from)
  upp = max(reference_coordinates$to)


  #change the solid and dashed lines for better separating chromosomes.
  p1 = ggplot2::ggplot(reference_coordinates) +
    CNAqc:::my_ggplot_theme(cex = cex) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = centromerStart,
        xend = centromerEnd,
        y = 0,
        yend = Inf
      ),
      size = .1,
      color = 'black',
      linetype = 8
    )

  p1 = p1 + ggplot2::geom_rect(
    data = reference_coordinates,
    ggplot2::aes(
      xmin = from,
      xmax = from,
      ymin = 0,
      ymax = Inf
    ),
    alpha = 1,
    colour = 'grey',
  )

  p1 = p1 +
    ggplot2::geom_hline(yintercept = 0,
                        size = 1,
                        colour = 'gainsboro') +
    ggplot2::geom_hline(
      yintercept = 1,
      size = .3,
      colour = 'black',
      linetype = 'dashed'
    ) +
    ggplot2::labs(x = "Chromosome",
                  y = "Major/ minor allele") +
    ggpubr::rotate_y_text() +
    # ggpubr::rotate_x_text() +
    # xlim(low, upp) +

    #set the chr names in the centromer positions.
    ggplot2::scale_x_continuous(
      breaks = c(0, reference_coordinates$centromerStart, upp),
      labels = c("", gsub(pattern = 'chr', replacement = '', reference_coordinates$chr), "")
    )

  return(p1)
}


add_segments_to_plot = function(segments, base_plot, cn)
{
  if (cn == 'absolute')
  {
    # Add one Major and minor lines, one on top of the other, red and blu
    M_seg = segments %>% dplyr::select(from, to, Major) %>% dplyr::rename(value = Major)
    m_seg = segments %>% dplyr::select(from, to, minor) %>% dplyr::rename(value = minor)

    base_plot = base_plot +
      ggplot2::geom_segment(
        data = M_seg %>% dplyr::mutate(Allele = "Major allele (clonal)"),
        ggplot2::aes(
          x = from,
          xend = to,
          y = value,
          yend = value,
          colour = Allele
        ),
        size = 1.5
      ) +
      ggplot2::geom_segment(
        data = m_seg %>% dplyr::mutate(Allele = "minor allele (clonal)"),
        ggplot2::aes(
          x = from,
          xend = to,
          y = value - 0.1,
          yend = value - 0.1,
          colour = Allele
        ),
        size = 1
      ) +
      ggplot2::scale_color_manual(values = c(`Major allele (clonal)` = 'red', `minor allele (clonal)` = 'steelblue')) +
      ggplot2::guides(color = ggplot2::guide_legend(''))

    # Some layout
    base_plot = base_plot +
      ggplot2::theme(
        legend.position = "bottom",
        legend.justification = "right",
        legend.margin = ggplot2::margin(0, 0, 0, 0)
      ) +
      ggplot2::labs(y = "Absolute allele counts")

  }

  if (cn == 'total')
  {
    base_plot = base_plot +
      ggplot2::geom_segment(
        data = segments %>% dplyr::select(from, to, total) %>% dplyr::mutate(Allele = "Segment ploidy"),
        ggplot2::aes(
          x = from,
          xend = to,
          y = total,
          yend = total
        ),
        size = 1.5,
        colour = 'black'
      )
  }

  return(base_plot)
}


plot_vaf = function(x, K){

  
  utils::globalVariables(c("..count.."))
  
  
  mutations <- x$mutations
  results <- x$results_timing
  reference_genome <- get_reference(x$reference_genome)
  purity = x$purity

  vfrom = reference_genome$from
  names(vfrom) = reference_genome$chr

  absolute_mutations <- mutations  %>%
    dplyr::mutate(from = .data$from + vfrom[.data$chr],
                  to = .data$to + vfrom[.data$chr])

  segments <- results$data$accepted_cna
  segments$from = lapply(segments$segment_name, function(s) {unlist(strsplit(s, "_"))[2]}) %>% unlist() %>% as.numeric()
  segments$to = lapply(segments$segment_name, function(s) {unlist(strsplit(s, "_"))[3]}) %>% unlist() %>% as.numeric()

  absolute_segments <- segments %>%
    dplyr::mutate(from = .data$from + vfrom[.data$chr],
                  to = .data$to + vfrom[.data$chr])

  summarized_results <- results$draws_and_summary[[K]]$summarized_results %>%
    dplyr::mutate(from = absolute_segments[.data$segment_id,]$from) %>%
    dplyr::mutate(to = absolute_segments[.data$segment_id,]$to) %>%
    dplyr::mutate(tau_mean = ifelse(.data$clock_mean  < 1, .data$clock_mean , 1)) %>%
    dplyr::mutate(tau_high = ifelse(.data$clock_high < 1, .data$clock_high, 1)) %>%
    dplyr::mutate(tau_low = .data$clock_low)



  accepted_mutations = data.frame()
  for (segment_idx in 1:nrow(summarized_results)) {
    segment <- summarized_results[segment_idx, ]
    segment_mutations <- absolute_mutations %>%
      dplyr::filter(.data$chr == segment$chr, .data$from > segment$from, .data$to < segment$to) %>%
      tidyr::drop_na(DP)
    # print(nrow(segment_mutations))
    cluster = summarized_results$tau_mean[segment_idx] # tau mean o clock mean?
    peaks = x$results_timing$data$input_data$peaks[[segment_idx]]
    # if (nrow(segment_mutations)> 40){
    segment_mutations <- segment_mutations %>% mutate(cluster = cluster)

    alpha = .05
    probs <- c(alpha/2, 1 - alpha/2)
    peaks = tickTack:::get_clonal_peaks(unique(segment_mutations$karyotype), purity)
    accepted_idx <- lapply(1:nrow(segment_mutations), function(i) {

      for (p in peaks) {
        quantiles <- stats::qbinom(probs, segment_mutations$DP[i], p)
        if ((segment_mutations$NV[i] >= quantiles[1]) && (segment_mutations$NV[i] <= quantiles[2])) {
          return(i)
        }
      }
    }) %>% unlist()
    segment_mutations <- segment_mutations[accepted_idx,]

    accepted_mutations <- bind_rows(accepted_mutations, segment_mutations)
    # }
  }
  as.factor(accepted_mutations$cluster)
  accepted_mutations <- accepted_mutations %>% mutate(VAF = NV/DP)

  accepted_mutations$cluster <- as.numeric(factor(accepted_mutations$cluster))

  k_colors = list(
    '2:0' = 'turquoise4',
    '2:1' = ggplot2::alpha('orange', .8),
    '2:2' = 'firebrick3'
  )
  my_palette <- c("#66a61e","#7570b3","#e7298a", "#1b9e77", "#d95f02")

  comb = expand_grid(unique(accepted_mutations$cluster), unique(accepted_mutations$karyotype))

  DENS = dplyr::tibble()
  PLOT_DATA = dplyr::tibble()
  for (i in 1:nrow(comb)) {
    cl = as.numeric(comb[i,1])
    kar = as.character(comb[i,2])

    df = accepted_mutations %>%
      dplyr::filter(cluster == cl, karyotype == kar) %>%
      dplyr::select(DP, NV, VAF) %>%
      dplyr::rename(successes = NV, trials = DP)

    if (nrow(df) > 0) {
      mobfit = mobster::mobster_fit(df, K = c(1,2), tail = FALSE, auto_setup = "FAST")

      binwidth = 0.01
      domain = seq(0, 1, binwidth)
      plot_data = mobster:::Clusters(mobfit$best)
      clusters = sort(unique(plot_data$cluster), na.last = TRUE)
      Beta_peaks = mobfit$best$Clusters %>%
        dplyr::filter(type == 'Mean', cluster != 'Tail')

      densities = suppressWarnings(mobster:::template_density(
        mobfit$best,
        x.axis = domain,
        binwidth = binwidth,
        reduce = TRUE
      ))

      DENS = dplyr::bind_rows(DENS, densities %>% dplyr::mutate(tickTack_cl = cl, karyotype = kar))
      PLOT_DATA = dplyr::bind_rows(PLOT_DATA, plot_data %>% dplyr::mutate(tickTack_cl = cl, karyotype = kar))
    }
  }

  COL = c("C1" = "#E41A1C", "C2" = "#377EB8")
  PLOT_DATA = PLOT_DATA %>%
    dplyr::group_by(karyotype, tickTack_cl) %>%
    dplyr::mutate(dens = VAF / sum(VAF))

  PLOT_DATA = PLOT_DATA %>% dplyr::mutate(cluster = ifelse(cluster == "C2", "P1", "P2"))
  DENS = DENS %>% dplyr::mutate(cluster = ifelse(cluster == "C2", "P1", "P2"))

  ggplot(PLOT_DATA, aes(VAF, fill = factor(tickTack_cl), y = ggplot2::after_stat(.data$..count.. / tapply(.data$..count.., .data$PANEL, sum)[.data$PANEL]))) +
    geom_histogram(alpha = 0.5, color = NA, position = "identity", binwidth = binwidth) +
    geom_line(data = DENS, aes(y = y, x = x, linetype=factor(cluster)), size = .5) +
    scale_color_manual(values = COL) +
    scale_fill_manual(values = my_palette) +
    theme_minimal() +
    scale_x_continuous(breaks = scales::pretty_breaks(n=3)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n=3)) +
    labs(x = "VAF", y = "Density", fill = "", col = "Cluster", linetype="Peak") +
    ggh4x::facet_nested(tickTack_cl~"CN"+karyotype) +
    # facet_grid(tickTack_cl ~ karyotype,
    #            labeller = labeller(var1 = label_both, var2 = label_both)) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2),
      strip.background = element_rect(fill = "gray80", color = "gray30")
    )
}