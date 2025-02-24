library(ggplot2)
library(plotly)
library(ggpubr)

#' Create annotations for plots based on subclonality counts
#' @param subclonality_counts Data frame containing counts of clonal/subclonal mutations
#' @return List of ggplot2 annotation objects
create_annotations <- function(subclonality_counts) {
  annotations <- list()
  
  if ("Clonal" %in% subclonality_counts$Subclonality) {
    n_clonal <- subclonality_counts$n[subclonality_counts$Subclonality == "Clonal"]
    if (length(n_clonal) > 0 && n_clonal > 0) {
      annotations <- c(
        annotations, 
        annotate("text", 
                 x = 0.85, 
                 y = -1.65, 
                 label = paste0("n clonal = ", n_clonal),
                 size = 3, 
                 color = "black", 
                 vjust = -1.5)
      )
    }
  }
  
  if ("Subclonal" %in% subclonality_counts$Subclonality) {
    n_subclonal <- subclonality_counts$n[subclonality_counts$Subclonality == "Subclonal"]
    if (length(n_subclonal) > 0 && n_subclonal > 0) {
      annotations <- c(
        annotations, 
        annotate("text", 
                 x = 0.85, 
                 y = -2.65, 
                 label = paste0("n subclonal = ", n_subclonal),
                 size = 3, 
                 color = "black", 
                 vjust = -1.5)
      )
    }
  }
  
  return(annotations)
}

#' Create VAF (Variant Allele Frequency) plot
#' @param merged_df Processed mutation data
#' @param sample_name Name of the sample
#' @param purity Purity value
#' @param read_depth List containing mean, min, and max read depth
#' @param smoothing_factor Smoothing factor for density plot
#' @param color_by Whether to color by "clonality" or "depth"
#' @param subclonality_counts Data frame with mutation counts
#' @param driver_mutations Optional data frame of driver mutations to highlight
#' @return plotly object
create_vaf_plot <- function(merged_df, sample_name, purity, read_depth, 
                            smoothing_factor, color_by = "clonality",
                            subclonality_counts, driver_mutations = NULL) {
  
  # Determine color aesthetic
  color_aesthetic <- if(color_by == "clonality") "Subclonality" else "ReadDepth"
  
  # Base plot
  plot <- ggplot(merged_df, aes(x = AF)) +
    geom_density(fill = "gray", alpha = .6, adjust = smoothing_factor) +
    geom_point(
      shape = 142,
      size = 5,
      alpha = 0.3,
      aes(
        y = -0.75,
        color = !!sym(color_aesthetic)
      )
    ) +
    xlim(0, 1) +
    ggtitle(paste(
      sample_name,
      " - Purity =", round(purity, 3),
      "      ", "Mean read depth =", round(read_depth$mean),
      "      ", "Range read depth =", read_depth$min, "-", read_depth$max
    )) +
    create_annotations(subclonality_counts) +
    ggpubr::theme_pubr() +
    xlab("VAF") + 
    ylab("Density") +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  
  # Add driver mutations if provided
  if (!is.null(driver_mutations)) {
    plot <- plot +
      geom_point(
        data = driver_mutations,
        aes(x = VAF, y = -0.75, fill = Subclonality),
        shape = 21,
        size = 8,
        color = "black",
        alpha = 0.7
      ) +
      geom_text(
        data = driver_mutations,
        aes(x = VAF, y = -1.25, label = gene),
        size = 3,
        angle = 45,
        hjust = 0
      )
  }
  
  # Apply color scheme
  if (color_by == "clonality") {
    plot <- plot +
      scale_color_manual(values = c("Clonal" = "#87CEEB", "Subclonal" = "#FFA07A")) +
      scale_fill_manual(values = c("Clonal" = "#87CEEB", "Subclonal" = "#FFA07A"))
  }
  
  # Convert to plotly
  ggplotly(plot)
}

#' Create CP (Cellular Prevalence) plot
#' @param merged_df Processed mutation data
#' @param sample_name Name of the sample
#' @param smoothing_factor Smoothing factor for density plot
#' @param color_by Whether to color by "clonality" or "depth"
#' @param subclonality_counts Data frame with mutation counts
#' @param df_subc Subclonal structure data
#' @param driver_mutations Optional data frame of driver mutations to highlight
#' @return plotly object
create_cp_plot <- function(merged_df, sample_name, smoothing_factor,
                           color_by = "clonality",
                           subclonality_counts, df_subc,
                           driver_mutations = NULL) {
  
  # Determine color aesthetic
  color_aesthetic <- if(color_by == "clonality") "Subclonality" else "ReadDepth"
  
  # Base plot
  plot <- ggplot(merged_df, aes(x = CP_unpenalized)) +
    geom_density(fill = "seagreen", alpha = .6, adjust = smoothing_factor) +
    geom_point(
      shape = 142,
      size = 5,
      alpha = 0.3,
      aes(
        y = -0.75,
        color = !!sym(color_aesthetic)
      )
    ) +
    xlim(0, 1.5) +
    ggtitle(paste(sample_name, " - Estimated mutation-specific CP")) +
    create_annotations(subclonality_counts) +
    ggpubr::theme_pubr() + 
    xlab("Estimated mutation-specific CP") + 
    ylab("Density") +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    geom_segment(
      data = df_subc, 
      aes(x = cellular_prevalence, 
          xend = cellular_prevalence, 
          y = 0, 
          yend = 5.7), 
      lty = "dashed"
    )
  
  # Add driver mutations if provided
  if (!is.null(driver_mutations)) {
    plot <- plot +
      geom_point(
        data = driver_mutations,
        aes(x = CP_unpenalized, y = -0.75, fill = Subclonality),
        shape = 21,
        size = 8,
        color = "black",
        alpha = 0.7
      ) +
      geom_text(
        data = driver_mutations,
        aes(x = CP_unpenalized, y = -1.25, label = gene),
        size = 3,
        angle = 45,
        hjust = 0
      )
  }
  
  # Apply color scheme
  if (color_by == "clonality") {
    plot <- plot +
      scale_color_manual(values = c("Clonal" = "#87CEEB", "Subclonal" = "#FFA07A")) +
      scale_fill_manual(values = c("Clonal" = "#87CEEB", "Subclonal" = "#FFA07A"))
  }
  
  # Convert to plotly
  ggplotly(plot)
}

#' Get VAF plot caption
#' @return HTML formatted caption for VAF plot
get_vaf_plot_caption <- function() {
  HTML("<p>
    - <b>Variant Allele Frequency (VAF)</b> indicates the proportion of sequencing reads that support a variant allele.<br>
    - <b>Purity</b> reflects the proportion of tumor cells in the sample. For instance, a purity of 0.9 means 90% of the cells are cancerous.<br>
    - <b>Mean read depth</b> is the average number of sequencing reads per SNV, influencing data reliability.<br>
    - <b>Density</b> estimates the probability density function of VAFs based on kernel density estimation (KDE). <br>
  </p>")
}

#' Get CP plot caption
#' @return HTML formatted caption for CP plot
get_cp_plot_caption <- function() {
  HTML("<p>
    - <b>Cellular Prevalence (CP)</b> the fraction of all cells (both tumor and admixed normal cells) from the sequenced tissue carrying a particular SNV.<br>
    - Annotations on the x-axis indicate whether mutations are <b>clonal</b> (present in all tumor cells) or <b>subclonal</b> (found in a subset).<br>
    - Counts of clonal and subclonal mutations are indicated as <b>n clonal</b> and <b>n subclonal</b> respectively.<br>
  </p>")
}

#' Configure plotly options for interactive plots
#' @param p plotly object
#' @return Configured plotly object
configure_plotly <- function(p) {
  p %>% config(
    displayModeBar = TRUE,
    scrollZoom = TRUE,
    toImageButtonOptions = list(
      format = "png",
      filename = "plot",
      height = 500,
      width = 700,
      scale = 2
    )
  )
}