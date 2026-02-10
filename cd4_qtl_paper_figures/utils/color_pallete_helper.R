#helper for color palletes

## plattes

ancestry_colors<-c("AGING COHORT"="#1b9e77",
                   "AFR"="#F2F1F0", 
                   "EUR"="#C0C9CC", 
                   "EAS"="#F9F6F0", 
                   "AMR"="#E8D4C0", 
                   "SAS"="#DAE7DD")


cohort_colors<-c("Female"="#81489c",
                 "Male"="#27bec9")

caQTL_category_cmap <- c(
  "in_caPeak" = "#b9dca4",
  "in_corr_Peak" = "#9dcdc1",
  "in_other_Peak" = "#099099",
  "no_Peak_overlap" = "#b1bcb5"
)
caQTL_category_cmap2 <- c(
  "in_caPeak" = "#b9dca4",
  "in_other_Peak" = "#9dcdc1",
  "no_Peak_overlap" = "#099099"
)

qtl_cmap<-c(
  "eQTLs"="#c7197c",
  "caQTLs"="#9ccb86"
)

variant_cmap2 <- c(
  "caQTL_only"  = "#9ccb86",
  "eQTL_only" = "#c7197c",
  "chromBPnet" = "#0071ba",
  "shared"     = "#f39c12"
)

variant_cmap3 <- c(
  "ChromBPNet_spec_only" = "#0071ba",
  "QTL_shared_any"     = "#f39c12"
)

variant_cmap4<- c(
  "ChromBPNet_spec" = "#0071ba",
  "QTL_shared"     = "#f39c12"
)
variant_cmap<-c(
  "eQTLs"="#c7197c",
  "caQTLs"="#9ccb86",
  "chromBPnet"="#0071ba"
)

qtl_cmap_coloc<-c(
  "just_eqtl"="#c7197c",
  "just_caqtl"="#9ccb86",
  "both" = "#f39c12")

qtl_cmap_coloc2<-c(
  "just_eqtl"="#c7197c",
  "just_caqtl"="#9ccb86",
  "both" = "#f39c12",
  "non_coloc"="gray90")

genotype_cmap<-c(
  "Homozygous Ref" = "darkgreen",
  "Heterozygous" = "gray40",
  "Homozygous Alt" = "darkviolet"
)


ld_colors <- c(
  "0.0–0.2" = "gray",
  "0.2–0.4" = "gray",
  "0.4–0.6" = "green",
  "0.6–0.8" = "orange",
  "0.8–1.0" = "red"
)

ld_colors2 <- c(
  "0.0–0.2" = "darkblue",
  "0.2–0.4" = "skyblue",
  "0.4–0.6" = "yellow",
  "0.6–0.8" = "orange",
  "0.8–1.0" = "red"
)


trait_labels <- c(
  "mono-count"   = "Monocyte Count",
  "lymph-count"  = "Lymphocyte Count",
  "eos-count"    = "Eosinophils Count",
  "IBD"          = "Inflammatory Bowel Disease",
  "AITD"         = "Autoimmune Thyroid Disease",
  "AST"          = "Asthma",
  "neutro-count" = "Neutrophil Count",
  "myeloid-count"= "Myeloid Count",
  "RA"           = "Rheumatoid Arthritis",
  "T1D"          = "Type 1 Diabetes",
  "SLE"          = "Systemic Lupus Erythematosus",
  "eos-pct"      = "Eosinophil %",
  "mono-pct"     = "Monocyte %",
  "lymph-pct"    = "Lymphocyte %",
  "CD"           = "Crohn's Disease",
  "neutro-pct"   = "Neutrophil %",
  "UC"           = "Ulcerative Colitis",
  "baso-count"   = "Basophil Count",
  "MS"           = "Multiple Sclerosis",
  "PSO"          = "Psoriasis",
  "AT"           = "Alergic Disease",
  "PDAST"        = "Pediatric Asthma",
  "VIT"          = "Vitiligo",
  "baso-pct"     = "Basophil %",
  "SS"           = "Systemic Sclerosis",
  "GD"           = "Graves’ Disease",
  "AD"           = "Contact Dermatitis",
  "PBC"          = "Primary Biliary Cirrhosis",
  "SD"           = "Sjögren's syndrome",
  "ECZ"          = "Eczema",
  "ADD"          = "Addison’s Disease",
  "MN"           = "Membranous Nephropathy",
  "MG"           = "Myasthenia Gravis",
  "HT"           = "Hashimoto’s Thyroiditis",
  "IL27"         = "IL27 Measurement",
  "JIA"          = "Juvenile Idiopathic Arthritis",
  "PV"           = "Psoriasis Vulgaris",
  "TNFR2"        = "TNFR2 Levels",
  "proIL16"      = "Pro-IL16 Measurement",
  "IL1"          = "IL1R antagonist Measurement",
  "SIGD"         = "Selective IgA Deficiency",
  "PSC"          = "Primary Sclerosing Cholangitis",
  "TNFR1"        = "TNFR1 Levels",
  "CCL4"         = "CCL4 Chemokine Measurement",
  "IL6"          = "IL6 Measurement",
  "LADA"="Latent Autoimmune Diabetes",
  "MCSF"="Macrophage Colony Stimulating Factor Measurement",
  "IL18"="IL18 Measurement",
  "CXCL1"="CXCL1 Measurement",
  "IGN"="IgA Nephropathy",
  "TNF"="TNF-related Apoptosis-inducing Ligand Measurement",
  "AS"="Ankylosing spondylitis"
  
  
)




