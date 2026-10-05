#turn off messages and warnings and make it so output isn't prefixed by anything,
#default is to put "##" in front of all output for some reason
#also set tidy to true so code is wrapped properly 
knitr::opts_chunk$set(message=FALSE, warning=FALSE, comment = "", cache = F, fig.align = "center")
options(width = 200)
# devtools::install_github("hoesler/rwantshue")
packagesUsed = c("GGally", "DT", "ggplot2", "stringr",
                 "tidyverse", "stringr", 'bookdown', 
                 'knitr', 'rmarkdown', "finalfit",
                 "modelsummary", "heatmaply", 
                 "fastmatch", "rwantshue", "lubridate", 
                 "plotly", "ggtree", "ComplexHeatmap", 
                 "dendextend", "magrittr", "ggpmisc", 
                 "ggthemes", "ggridges", "here")
library(rwantshue)
library(ggthemes)
library(here)
suppressMessages(lapply(packagesUsed, require, character.only = TRUE))
# devtools::install_github("hoesler/rwantshue")
myFormula= x~y
# commonly used functions  
`%!in%` <- Negate(`%in%`)
is.notna <-function(x){
  return(!is.na(x))
}
# create a color scheme object
scheme <- iwanthue() # get the iwanthue singleton
scheme <- iwanthue(seed = 42, force_init = TRUE) # recreate with a seed

# generate a new color palette (vector of hex values) with presets...
scheme$hex(10)
scheme$hex(10, color_space = hcl_presets$fluo)

# ggplot themes 
transparentBackground = theme(
  panel.background = element_rect(fill='transparent'), #transparent panel bg
  plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
  panel.grid.major = element_blank(), #remove major gridlines
  panel.grid.minor = element_blank(), #remove minor gridlines
  legend.background = element_rect(color = NA, fill='transparent'), #transparent legend bg
  legend.box.background = element_rect(color = NA, fill='transparent') #transparent legend panel
)
sofonias_theme_noTransparentBackground = theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank() )+
  theme(axis.line.x = element_line(color="black", linewidth = 0.3),axis.line.y =
          element_line(color="black", linewidth = 0.3))+
  theme(text=element_text(size=12, family="Helvetica"))+
  theme(axis.text.y = element_text(size=12))+
  theme(axis.text.x = element_text(size=12)) +
  theme(legend.position = "bottom") + 
  theme(plot.title = element_text(hjust = 0.5))
sofonias_theme = theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank() )+
  theme(axis.line.x = element_line(color="black", linewidth = 0.3),axis.line.y =
          element_line(color="black", linewidth = 0.3))+
  theme(text=element_text(size=12, family="Helvetica"))+
  theme(axis.text.y = element_text(size=12))+
  theme(axis.text.x = element_text(size=12)) +
  theme(legend.position = "bottom") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  transparentBackground

sofonias_theme_xRotate_noTransparentBackground = sofonias_theme + 
  theme(axis.text.x = element_text(size=12, angle = -90, vjust = 0.5, hjust = 0)) 
sofonias_theme_xRotate = sofonias_theme + 
  theme(axis.text.x = element_text(size=12, angle = -90, vjust = 0.5, hjust = 0)) + 
  transparentBackground


# color schemes

colorPalette_08 = c("#2271B2","#3DB7E9","#F748A5","#359B73","#D55E00","#E69F00","#F0E442","#000000")
colorPalette_12 = c("#9F0162","#009F81","#FF5AAF","#00FCCF","#8400CD","#008DF9","#00C2F9","#FFB2FD","#A40122","#E20134","#FF6E3A","#FFC33B")
colorPalette_15 = c("#68023F","#008169","#EF0096","#00DCB5","#FFCFE2","#003C86","#9400E6","#009FFA","#FF71FD","#7CFFFA","#6A0213","#008607","#F60239","#00E307","#FFDC3D")

mycol1 = rgb(85, 138, 221, maxColorValue = 255)
mycol2 = rgb(255, 192, 78, maxColorValue = 255)
mycol3 = rgb(255, 121, 177, maxColorValue = 255)
mycol4 = rgb(221, 221, 221, maxColorValue = 255)
mycol5 = rgb(100, 89, 89, maxColorValue = 255)
mycol6 = rgb(0, 189, 189, maxColorValue = 255)

parasiteDensityColors = c("10" = mycol2, "100" = mycol3, "1000" = mycol6, "10000" = mycol1)

heome1_targetsToRemove = c("t2","t21","t40","t41","t56","t64","t72","t76","t93")


templateToDensity = tibble(s_startingTemplate = c("ST0006", "ST0060", "ST0600", "ST6000"), 
                           s_parasiteDensity = c(10, 100, 1000, 10000))


palettes <- ggthemes_data[["tableau"]][["color-palettes"]][["regular"]]

create_dt <- function(x) {
  DT::datatable(
    x,
    extensions = 'Buttons',
    options = list(
      dom = 'Blfrtip',
      buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
      lengthMenu = list(c(10, 25, 50,-1),
                        c(10, 25, 50, "All"))
    ), 
    filter = "top"
  )
}
programColors = c("SeekDeep" = "#2271B2", "dada2" = "#D55E00")



createColorListFromDf <- function(df, colorPalette = colorPalette_12, iwanthudSeed = rnorm(1) * 100){
  colorList = list()
  for (dfColname in colnames(df)) {
    levels = sort(unique(df[[dfColname]]))
    scheme <- iwanthue(seed = iwanthudSeed, force_init = TRUE)
    if (length(levels) <= length(colorPalette_12)) {
      levelsCols = colorPalette_12[1:length(levels)]
    } else{
      levelsCols = scheme$hex(length(levels))
    }
    names(levelsCols) = levels
    colorList[[dfColname]] = levelsCols
  }
  return(colorList)
}

createDownloadLink <-function(fnp, linkName = ""){
  relToHere = sub(here(),"",normalizePath(getwd()), fixed = T)
  depth = stringr::str_count(relToHere, "/")
  prefix = ""
  if(depth > 0){
    prefix = paste0(paste0(rep("..",depth), collapse = "/"), "/")
  }
  if ("" == linkName) {
    linkName =  basename(fnp)
  }
  return(paste0("[", linkName, "]", "(", paste0(prefix, gsub(paste0(here(),"/"),"",normalizePath(fnp), fixed = T)), "){.downloadLink .btn .btn-info}" ) )
}

createImgLink <-function(fnp, linkName = ""){
  relToHere = sub(here(), "", normalizePath(getwd()), fixed = T)
  depth = stringr::str_count(relToHere, "/")
  prefix = ""
  if(depth > 0){
    prefix = paste0(paste0(rep("..",depth), collapse = "/"), "/")
  }
  if ("" == linkName) {
    linkName =  basename(fnp)
  }
  return(paste0("![", linkName, "]", "(", paste0(prefix, gsub(paste0(here(),"/"),"",normalizePath(fnp), fixed = T)), ")" ) )
}

create_dt <- function(x){
  DT::datatable(x,
                extensions = 'Buttons',
                options = list(dom = 'Blfrtip',
                               buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                               lengthMenu = list(c(10,25,50,-1),
                                                 c(10,25,50,"All"))), 
                filter = "top")
}



create_tabsetOfHtmlWidgets <- function(htmlObjectsList) {
  zz <- textConnection("foo", "w")
  sink(zz)
  cat("::: {.panel-tabset}\n")
  iwalk(htmlObjectsList, ~ {
    cat('## ', .y, '\n\n')
    tempList = list()
    tempList[["item"]] = .x
    print(htmltools::tagList(tempList))
    cat('\n\n')
  })
  cat(":::\n")
  sink()
  close(zz)
  paste0(foo, collapse = "\n")
}

create_tabsetOfGgplotObjects <- function(ggplotObjectsList) {
  # zz <- textConnection("foo", "w")
  # sink(zz)
  cat("\n\n::: {.panel-tabset}\n\n")
  iwalk(ggplotObjectsList, ~ {
    cat('## ', .y, '\n\n')
    print(.x)
    # tempList = list()
    # tempList[["item"]] = .x
    # grid::grid.draw(ggplot2::ggplotGrob(.x))
    cat('\n\n')
  })
  cat(":::\n")
  # sink()
  # close(zz)
  # paste0(foo, collapse = "\n")
}

writeFasta<-function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"name"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"seq"]))
  }
  if(stringr::str_ends(filename, ".gz")){
    fileConn<-gzfile(filename)
    writeLines(fastaLines, fileConn)
    close(fileConn)
  }else{
    fileConn<-file(filename)
    writeLines(fastaLines, fileConn)
    close(fileConn) 
  }
}


genShapesExample <- function() {
  df_shapes <- data.frame(shape = 0:24)
  return(  ggplot(df_shapes, aes(0, 0, shape = shape)) +
             geom_point(aes(shape = shape), size = 5, fill = 'red') +
             scale_shape_identity() +
             facet_wrap( ~ shape) +
             theme_void())
}

genColorsForValues <- function(values) {
  uniqValues = unique(values)
  outColors = c()
  if (length(uniqValues) <= 8) {
    outColors = colorPalette_08[1:length(uniqValues)]
  } else if (length(uniqValues) <= 12) {
    outColors = colorPalette_12[1:length(uniqValues)]
  } else if (length(uniqValues) <= 15) {
    outColors = colorPalette_15[1:length(uniqValues)]
  } else{
    outColors = scheme$hex(length(uniqValues))
  }
  names(outColors) = uniqValues
  return (outColors)
}


#' Check sharing and unqiue values between two vectors
#'
#' @param vectorA the first vector
#' @param vectorB the second vector
#'
#' @return a list with 4 vectors, "only_in_vectorA" unique to vectorA, "only_in_vectorB" unique to vectorB, "shared_samples" shared between both vectorA and vectorB, "all" all values by combinng vectorA and vectorB 
set_decompose <- function(vectorA, vectorB){
  ret = list()
  # Find unique and shared samples
  ret[["only_in_vectorA"]] <- setdiff(vectorA, vectorB)  # Samples only in vectorA
  ret[["only_in_vectorB"]] <- setdiff(vectorB, vectorA)  # Samples only in vectorB
  ret[["shared_samples"]] <- intersect(vectorA, vectorB) # Samples shared between vectorA and vectorB
  ret[["all"]] <- union(vectorA, vectorB) # All samples between vectorA and vectorB
  return(ret)
}

# define once
quick_summary <- function(x, probs = c(0.25, 0.5, 0.75)) {
  q <- quantile(x, probs = probs, na.rm = TRUE)
  tibble(
    n = length(x[!is.na(x)]), 
    sd = sd(x, na.rm = TRUE),
    mean = mean(x, na.rm = TRUE),
    min  = min(x, na.rm = TRUE),
    max  = max(x, na.rm = TRUE),
    !!!setNames(as.list(q), paste0("q", probs*100))
  )
}

#' Add various contingency table metrics, input must have TP, FP, TN, FN columns defined already
#'
#' @param input A data table with at least TP, FP, TN, FN
#' @param summarize_first Whether or not to sum up the columns first, will lose any non-grouping columns though
#'
#' @returns the input table plus various performance metrics 
#' @export
add_contigency_metrics <- function(input, summarize_first = F){
  output = input
  if(summarize_first){
    output = output  %>%
      summarise(TP = sum(TP), 
                FN = sum(FN), 
                FP = sum(FP), 
                TN = sum(TN))
  }
  output = output %>%
    mutate(
      sensitivity           = if_else(TP + FN > 0, TP / (TP + FN), NA_real_),
      recall                = sensitivity,
      specificity           = if_else(TN + FP > 0, TN / (TN + FP), NA_real_),
      ppv                   = if_else(TP + FP > 0, TP / (TP + FP), NA_real_),
      precision             = ppv,
      npv                   = if_else(TN + FN > 0, TN / (TN + FN), NA_real_),
      accuracy              = (TP + TN) / (TP + TN + FP + FN),
      prevalence            = (TP + FN) / (TP + TN + FP + FN),
      overlap_coefficient   = ifelse(TP + FP + FN > 0, TP / (TP + FP + FN), NA_real_),
      tpr                   = sensitivity,
      fpr                   = if_else(FP + TN > 0, FP / (FP + TN), NA_real_),
      fdr                   = if_else(TP + FP > 0, FP / (TP + FP), NA_real_),
      f1_score              = if_else((2 * TP + FP + FN) > 0, (2 * TP) / (2 * TP + FP + FN), NA_real_),
      balanced_accuracy     = (sensitivity + specificity) / 2
    )
  return(output)
}


#' Add various contingency table metrics but only if they don't depend on TN, input must have TP, FP, FN columns defined already
#'
#' @param input A data table with at least TP, FP, FN
#' @param summarize_first Whether or not to sum up the columns first, will lose any non-grouping columns though
#' @returns the input table plus various performance metrics 
#' @export
add_contigency_metrics_no_TN <- function(input, summarize_first = F){
  output = input
  if(summarize_first){
    output = output  %>%
      summarise(TP = sum(TP), 
                FN = sum(FN), 
                FP = sum(FP))
  }
  output = output %>%
    mutate(
      sensitivity           = if_else(TP + FN > 0, TP / (TP + FN), NA_real_),
      recall                = sensitivity,
      ppv                   = if_else(TP + FP > 0, TP / (TP + FP), NA_real_),
      precision             = ppv,
      overlap_coefficient   = ifelse(TP + FP + FN > 0, TP / (TP + FP + FN), NA_real_),
      tpr                   = sensitivity,
      fdr                   = if_else(TP + FP > 0, FP / (TP + FP), NA_real_),
      f1_score              = if_else((2 * TP + FP + FN) > 0, (2 * TP) / (2 * TP + FP + FN), NA_real_)
    )
  return(output)
}



