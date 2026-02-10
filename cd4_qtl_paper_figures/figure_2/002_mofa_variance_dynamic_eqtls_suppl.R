library(ggplot2)
library(MultiAssayExperiment)
library(MOFA2)
library(ggplot2)
library(MOFAdata)
library(cowplot)
options(bitmapType="cairo")



model34 <- load_model("/gchm/cd4_CellRegMap/001_preprocessing/results_01312025/mofa_factors/mofa_trained/cd4_aging_filt_sce_mofa_expectations_res34.hdf5")
model4 <- load_model("/gchm/cd4_CellRegMap/001_preprocessing/results_01312025/mofa_factors/mofa_trained/cd4_aging_filt_sce_mofa_expectations_res4.hdf5")
model45 <- load_model("/gchm/cd4_CellRegMap/001_preprocessing/results_01312025/mofa_factors/mofa_trained/cd4_aging_filt_sce_mofa_expectations_res45.hdf5")
model5 <- load_model("/gchm/cd4_CellRegMap/001_preprocessing/results_01312025/mofa_factors/mofa_trained/cd4_aging_filt_sce_mofa_expectations_res5.hdf5")

head(model34@cache$variance_explained$r2_per_factor) # group 1
head(model4@cache$variance_explained$r2_per_factor) # group 1
head(model45@cache$variance_explained$r2_per_factor) # group 1
head(model5@cache$variance_explained$r2_per_factor) # group 1

######### MOFA heatmap########
plot_variance_explained_composite <- function(models,
                                              model_names = NULL,
                                              factors = "all",
                                              min_r2 = 0,
                                              max_r2 = NULL,
                                              legend = TRUE,
                                              use_cache = TRUE,
                                              ...) {
  stopifnot(is.list(models))
  
  # If names not given, try from list; otherwise make generic ones
  if (is.null(model_names)) {
    model_names <- names(models)
    if (is.null(model_names) || any(model_names == "")) {
      model_names <- paste0("model", seq_along(models))
    }
  }
  stopifnot(length(model_names) == length(models))
  
  all_df <- Map(function(object, mname) {
    # 1) Get variance explained list (reuse MOFA2 logic)
    if (use_cache && .hasSlot(object, "cache") &&
        ("variance_explained" %in% names(object@cache))) {
      r2_list <- object@cache$variance_explained
    } else {
      r2_list <- calculate_variance_explained(object, factors = factors, ...)
    }
    
    r2_mk <- r2_list$r2_per_factor   # list over groups: each is factors x views
    
    # 2) Decide which factors to keep (same logic as original)
    if ((length(factors) == 1) && (factors[1] == "all")) {
      f_used <- factors_names(object)
    } else if (is.numeric(factors)) {
      f_used <- factors_names(object)[factors]
    } else {
      stopifnot(all(factors %in% factors_names(object)))
      f_used <- factors
    }
    
    # 3) Sum variance across views *and* groups for each factor:
    #    Sum all group matrices, then sum across columns (views)
    #    -> one value per factor
    mat_sum_groups <- Reduce("+", r2_mk)  # factors x views
    vec_factor_total <- rowSums(mat_sum_groups) # sum over views
    
    # subset to chosen factors and build df
    vec_factor_total <- vec_factor_total[f_used]
    
    data.frame(
      factor = f_used,
      model  = mname,
      value  = as.numeric(vec_factor_total),
      stringsAsFactors = FALSE
    )
  }, object = models, mname = model_names)
  
  all_df <- do.call(rbind, all_df)
  
  # 4) Tidy factor/model ordering
  all_df$factor <- factor(all_df$factor,
                          levels = rev(unique(all_df$factor)))  # top = Factor1
  all_df$model  <- factor(all_df$model, levels = model_names)
  
  # 5) Apply min/max clipping similarly to original
  if (!is.null(min_r2)) {
    all_df$value[all_df$value < min_r2] <- 0.001
    min_r2 <- 0
  } else {
    min_r2 <- min(all_df$value, na.rm = TRUE)
  }
  
  if (!is.null(max_r2)) {
    all_df$value[all_df$value > max_r2] <- max_r2
  } else {
    max_r2 <- max(all_df$value, na.rm = TRUE)
  }
  
  # 6) Plot heatmap: rows = factors, columns = models
  p <- ggplot(all_df, aes(x = model, y = factor, fill = value)) +
    geom_tile(color = "black") +
    scale_fill_gradientn(
      colors = c("gray97", "darkblue"),
      limits = c(min_r2, max_r2),
      guide = "colorbar"
    ) +
    guides(fill = guide_colorbar("Var. (%)")) +
    labs(x = "", y = "") +
    theme(
      axis.text.x = element_text(size = rel(1),  angle = 45, hjust = 1,
                                 color = "black"),
      axis.text.y = element_text(size = rel(1.1), color = "black"),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.background = element_blank(),
      axis.line = element_blank(),
      legend.title = element_text(color = "black"),
      legend.text  = element_text(color = "black")
    )
  
  if (isFALSE(legend)) {
    p <- p + theme(legend.position = "none")
  }
  
  return(p)
}
p_comp <- plot_variance_explained_composite(
  models = list(
    model34 = model34,
    model4  = model4,
    model45 = model45,
    model5  = model5
  ),
  factors = "all"   # or e.g. factors = 1:10
)

p_comp

pdf("~/cd4_qtl_paper_figures/supplements/plots/mofa_factors_pcs_all.pdf", width = 6, height = 4)
p_comp
dev.off()

