library(ggplot2)
library(ggseqlogo)
library(rtracklayer)
library(patchwork)
library(readr)
library(magrittr)
library(dplyr)

#####################
####################
### USAGE ##########
# #paths 
# prediction="~/cd4_chrombpnet/chrombpnet_model_b7/prediction_scores_bw/averaged/cd4_tcells_chrombpnet_nobias.bw"
# contribution="~/cd4_chrombpnet/chrombpnet_model_b7/contribution_scores_bw/averaged/cd4_tcells.counts_scores.bw"
# contribution2="~/cd4_chrombpnet/chrombpnet_model_b7/contribution_scores_bw/averaged/cd4_tcells.counts_scores.bw"
# 
# # Create a tibble
# config <- tibble(
#   sample=c("predicted", "contribution_bw", "contribution_dynaseq"),
#   path = c(prediction, contribution, contribution2),
#   normalize = c("TRUE", "TRUE", "TRUE"),
#   normalization_region = c('None', 'None', 'None'),
#   type = c('single', 'single', 'dynseq'),
#   Min = c('None', -0.03, -0.03),
#   Max =c('None', 0.03, 0.03),
#   Scale=c('3.26367792781148', '1.49254637494212', 'None'),
#   Color=c('1B9E77', '1B9E77', 'None'),
#   
# )
## Plot region
# range="chr2:135004329-135004583"
# cellline="cd4"
# REGION = GRanges(range)
# options(repr.plot.width = 40, repr.plot.height = nrow(config))
# p = generate_plots(config, REGION, genome)
# ggsave("plot_test.png", p, height = 4, width = 16)
# #################
################
#################

## Load Data
#Set the reference genome
library(BSgenome.Hsapiens.UCSC.hg38)
genome = BSgenome.Hsapiens.UCSC.hg38

# one bigwig at a time, any number of peaks
get_matrix_from_bigwig <- function(bigwig_path, peak_set) {
  # ensure peak set has fixed width
  stopifnot(length(unique(width(peak_set)))==1)
  
  as.matrix(import(bigwig_path, 
                   which=peak_set, as="NumericList"))
}

# one bigwig, one peak (GRanges object)
get_importance_from_bigwig <- function(bigwig_path, peak, genome) {
  stopifnot(length(peak)==1)
  
  # get DNA sequence
  sequence = genome[[as.vector(seqnames(peak))]][start(peak):end(peak)]
  
  m = matrix(0, length(sequence), 4)
  colnames(m) = c("A", "C", "G", "T")
  m[cbind(seq(length(sequence)), as.vector(matrix(sequence)))] = get_matrix_from_bigwig(bigwig_path, peak)
  
  t(m)
}

GRangesFromDataFrames<-function(dataframe){with(dataframe,GRanges(IRanges(start=start,end=end),seqnames=seqnames,strand=strand,dataframe%>%dplyr::select(-strand,-start,-end,-seqnames)))}

calculate_total_value_and_length <- function(bw_path,region){
  total_value = as.integer(sum(as.vector(get_matrix_from_bigwig(bw_path, region))))
  length = length(as.vector(get_matrix_from_bigwig(bw_path, region)))
  df = data.frame(total_value=total_value,length=length)
  df
}

calculate_normalization_value <- function(bw_path,peaks_path,mc.cores=40){
  peaks_df = read_tsv(peaks_path,
                      col_names=c("seqnames","start","end","name","score","strand","p","q","x","summit"),
                      show_col_types = FALSE)
  peaks_df['strand']<-'*'
  valid_chrs = c('chr1',"chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
                 "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20",
                 "chr21")
  peaks_df[peaks_df[["seqnames"]] %in% valid_chrs,]
  peaks_gr = GRangesFromDataFrames(peaks_df)
  
  peaks_gr = keepSeqlevels(peaks_gr, c('chr1',"chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
                                       "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20",
                                       "chr21"), pruning.mode="coarse")
  
  print(bw_path)
  1:length(peaks_gr) %>% mclapply(function(i){
    calculate_total_value_and_length(bw_path,peaks_gr[i])
  },mc.cores=mc.cores
  ) %>% bind_rows() -> normalization_df
  
  sum(normalization_df["total_value"])/sum(normalization_df["length"])
}

## Plot types
TITLE_SZ = 16

### Importance 
plot_seq <- function(m, ymin=NULL, ymax=NULL, clip=F) {
  mat = m
  
  # Cap values to ymax if clipping is enabled
  if (clip == TRUE) {
    if (!is.null(ymax)) {
      mat[mat > ymax] <- ymax  # Corrected the condition
    }
  }
  
  p = ggseqlogo(mat, method='custom', seq_type='dna') 
  
  p = (p + 
         theme_classic() + 
         coord_cartesian(ylim=c(ymin, ymax)) + 
         expand_limits(x=0, y=0) +
         theme(axis.ticks.x = element_blank(), 
               axis.text.x = element_blank(),
               axis.title.x = element_blank(),
               axis.line.x = element_blank()))
  
  p = p + theme(axis.ticks.y = element_blank(),
                axis.text.y = element_blank())
  
  p
}
plot_importance <- function(bigwig_path, region, genome, ylabel, ymin=NULL, ymax=NULL, clip=F) {
  plot_seq(get_importance_from_bigwig(bigwig_path, region, genome), ymin=ymin, ymax=ymax, clip=clip) +
    ylab(gsub("\\n", "\n", ylabel)) +
    theme(axis.title.y = element_text(angle=0, size=TITLE_SZ, hjust=0.5, vjust=0.5))
}
### Single Track
plot_single_vals <- function(v, ymin=NULL, ymax=NULL, 
                             clip=F, fill=T, 
                             col='#37ada2',x_width) {
  vals = v
  
  #     rownames(mat) = c("A", "C", "G", "T")
  
  # cap to upper and lower limits
  if (clip==T) {
    if (!is.null(ymin)) {
      vals[valsymax] = ymax
    }
  }
  
  p = ggplot(data.frame(x=seq(length(vals)), y=vals), aes(x=x,y=y)) 
  
  if (fill==T) {
    p = p + geom_area(fill=col)
  }
  
  else {
    p = p + geom_line(col=col)
  }
  
  p = (p + 
         theme_classic() + 
         coord_cartesian(ylim=c(0, ceiling(max(vals)))) + 
         expand_limits(x=0, y=0) +
         #scale_y_continuous(breaks=pretty(c(0, max(vals)),n=2,min.n=1)) +
         scale_y_continuous(breaks=c(0, ceiling(max(vals)))) +
         theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),axis.title.x = element_blank(),axis.line.x = element_blank()))
  
  p = p + theme(axis.text.y = element_text(size=14))
  
  p
}

plot_single <- function(bigwig_path, region, 
                        ylabel, ymin=NULL, ymax=NULL, 
                        fill=T, col='#37ada2', clip=F,
                        normalize = FALSE,normalization_value=None, normalization_region,normalization_file) {
  v1 = as.vector(get_matrix_from_bigwig(bigwig_path, region))
  if (normalize){
    #normalization_value = calculate_normalization_value(normalization_file,normalization_region)
    print(max(v1))
    print(normalization_value)
    v1=(v1/normalization_value)
    print(max(v1))
  }
  plot_single_vals(v1, ymin=ymin, ymax=ymax, clip=clip, fill=fill, col=col, x_width=width(region)) +
    ylab(gsub("\\n", "\n", ylabel)) +
    theme(axis.title.y = element_text(angle=0, size=TITLE_SZ, hjust=0.5,  vjust=0.5))
}

### scale
plot_scale_vals <- function(v, col='#a9d1ac') {
  vals = v
  
  p = ggplot(data.frame(x=seq(length(vals)), y=0), aes(x=x,y=y)) 
  
  
  p = p + geom_line(col=col)
  
  p = (p + 
         theme_classic() + 
         expand_limits(x=0, y=0) +
         scale_x_continuous(breaks=pretty(c(0, length(vals)/10),n=2,min.n=1)) +
         theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(),axis.title = element_blank(),axis.line.y = element_blank())
  ) 
  
  p
}

plot_scale <- function(bigwig_path, region, ylabel, ymin=NULL, ymax=NULL, fill=T, col='#37ada2', clip=F) {
  plot_scale_vals(as.vector(get_matrix_from_bigwig(bigwig_path, region)), col=col) +
    ylab(ylabel) +
    theme(axis.title.y = element_text(angle=0, size=TITLE_SZ, hjust=0.5, vjust=0.5))
}


### Stranded Track 
plot_stranded_vals <- function(v1, v2, ymin=NULL, 
                               ymax=NULL, clip=F, 
                               col1='blue', col2='orange') {
  vals1 = v1
  vals2 = v2
  
  #     rownames(mat) = c("A", "C", "G", "T")
  
  # cap to upper and lower limits
  if (clip==T) {
    if (!is.null(ymin)) {
      vals1[vals1ymax] = ymax
      vals2[vals2>ymax] = ymax
    }
  }
  
  p = ggplot(data.frame(x=seq(length(vals1)), y1=vals1, y2=vals2)) +
    geom_area(aes(x=x, y=y1), fill=col1, alpha = 0.8) +
    geom_area(aes(x=x, y=y2), fill=col2, alpha = 0.8) 
  
  p = (p + 
         theme_classic() + 
         coord_cartesian(ylim=c(ymin, ymax)) + 
         expand_limits(x=0, y=0) +
         scale_y_continuous(breaks=pretty(c(0, max(max(vals1),max(vals2))),n=2,min.n=1)) +
         theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),axis.title.x = element_blank(),axis.line.x = element_blank()))
  
  p
  
}

plot_stranded <- function(bigwig_prefix, region, ylabel, ymin=NULL, ymax=NULL,
                          clip=F, col1='blue', col2='orange',
                          normalize = FALSE,normalization_region,normalization_file_prefix) {
  # plus and minus should be bigwig_prefix + "plus.bw" and bigwig_prefix + "minus.bw"
  
  v1 = as.vector(get_matrix_from_bigwig(paste(bigwig_prefix, "pos.bw", sep='') , region))
  v2 = as.vector(get_matrix_from_bigwig(paste(bigwig_prefix, "neg.bw", sep='') , region))
  if (normalize){
    normalization_value = calculate_normalization_value(paste(normalization_file_prefix, "pos.bw", sep=''),normalization_region)
    v1=(v1/normalization_value)
  } 
  if (normalize){
    normalization_value = calculate_normalization_value(paste(normalization_file_prefix, "neg.bw", sep=''),normalization_region)
    v2=(v2/normalization_value)
  } 
  plot_stranded_vals(v1, v2, ymin=ymin, ymax=ymax, clip=clip, col1=col1, col2=col2) + 
    ylab(ylabel) + 
    theme(axis.title.y = element_text(angle=0, size=TITLE_SZ, hjust=0))
}

### Generate plots from config
#my_list <- list("#E00FEE", "#E00FEE", "#7B241C", "red", "#4A235A", "#4A235A", "#239B56", "blue")
#my_list <- list("#4A235A", "#4A235A", "#239B56", "blue")
#my_list <- list("#E00FEE", "#E00FEE", "null", "#458B00", "null", "#0000FF", "null", "#A52A2A", "null", "#A52A2A", "null", "#A52A2A", "null")

generate_plots <- function(config, region, genome) {
  plots = list()
  
  # configure for different relative heights for single, stranded and dynseq tracks
  REL_HEIGHTS = c(single=1, stranded=1, dynseq=1, scale =0.5)
  REL_widths = c(single=1, stranded=1, dynseq=1, scale =1)
  
  heights = c()
  
  for (i in seq(nrow(config))) {
    print(config$sample[i])
    if (config$type[i] == "dynseq") {
      plots[[i]] = plot_importance(config$path[i], region, genome, config$sample[i],
                                   ymin=as.double(config$Min[i]), 
                                   ymax=as.double(config$Max[i]), clip=T)            
    }
    
    else if (config$type[i] == "single") {
      #print(i)
      #print(my_list[i])
      #print(paste("#",config$Color[i],sep=""))
      plots[[i]] = plot_single(config$path[i], region, config$sample[i],
                               col=paste("#",config$Color[i],sep=""),
                               normalization_value=as.double(configScale[i],
                                                             normalize=config$normalize[i],
                                                             normalization_region=config$normalizationregion[i],
                                                             normalizationfile=config$normalization_file[i]
                               ))
    }
    else if (config$type[i] == "stranded") {
      plots[[i]] = plot_stranded(config$path[i], region, config$sample[i],
                                 normalize=confignormalize[i], normalization_region=config$normalization_region[i],
                                 normalization_file=config$normalization_file[i])
    }
    else if (config$type[i] == "scale") {
      plots[[i]] = plot_scale(config$path[i], region, config$sample[i])
    }
    
    heights = c(heights, REL_HEIGHTS[config$type[i]])        
    
  }
  
  main = wrap_plots(plots, heights=heights,  guides = "collect")  + 
    theme(plot.margin = margin(5, 5, 5, 5))
  main
}

