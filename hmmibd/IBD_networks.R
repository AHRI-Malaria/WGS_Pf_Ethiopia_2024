## script adapted from https://github.com/vivaxgen

library(igraph)
library(unikn)

setwd("~/Google Drive/My Drive/Research_Fellow/Popgen/ibd_networks")

mypch_shape <- function(pch_val) {
  function(coords, v = NULL, params) {
    vertex.color <- params("vertex", "color")
    if (length(vertex.color) != 1 && !is.null(v)) {
      vertex.color <- vertex.color[v]
    }
    vertex.size <- 1 / 200 * params("vertex", "size")
    if (length(vertex.size) != 1 && !is.null(v)) {
      vertex.size <- vertex.size[v]
    }
    points(coords[, 1], coords[, 2], pch = pch_val,
           bg = vertex.color, col = "black",
           cex = vertex.size * 10, lwd = 0.5)
  }
}

add_shape("triangle",
  clip = shapes("circle")$clip,
  plot = mypch_shape(24)
)

add_shape("diamond",
  clip = shapes("circle")$clip,
  plot = mypch_shape(23)
)


.reorder.metadata <-
  function(IBD, metadata, sample.col = "Sample") {
    samples <- unique(c(IBD[, "sample1"], IBD[, "sample2"]))
    sample.order <- match(metadata[, sample.col], samples)
    actual.length <- length(na.omit(sample.order))
    
    metadata <- metadata[order(sample.order), ]
    metadata[1:actual.length, ]
  }


.get.labels <-
  function(metadata,
           label.col,
           label.palette = NULL) {
    labels <- metadata[, label.col]
    labels[is.na(labels)] <- NaN
    
    labels <- factor(labels)
    
    if (is.null(label.palette)) {
      labels
    } else {
      factor(labels, levels = names(label.palette))
    }
  }


.generate.label.palette <- function(labels, palette = "custom") {
  label.names <- levels(labels)
  num.labels <- length(label.names)
  
  if (palette == "custom") {
    custom_palette <- c(
      "#A6CEE3", 
      "#238b8d", 
      "#a72e61", 
      "#FDBF6F", 
      "#CAB2D6", 
      "#FFFF99", 
      "#1F78B4", 
      "#33A02C", 
      "#E31A1C", 
      "#FF7F00", 
      "#6A3D9A", 
      "#B15928",
      "#F4A582", 
      "#92C5DE", 
      "#D9F0D3"  
    )
    
    if (num.labels > length(custom_palette)) {
      label.palette <- rep(custom_palette, length.out = num.labels)
    } else {
      label.palette <- custom_palette[1:num.labels]
    }
  } else if (palette == "unikn") {
    label.palette <- unikn::usecol(pal = pal_unikn, n = num.labels)
  } else {
    label.palette <- palette.colors(num.labels, palette = palette)
  }
  
  names(label.palette) <- label.names
  label.palette
}



.generate.label.colours <-
  function(labels, label.palette)
    label.palette[labels]


.create.edgelist <-
  function(IBD)
    IBD[, c("sample1", "sample2", "fract_sites_IBD")]


.create.vertices <-
  function(IBD,
           metadata = NULL,
           sample.col = "Sample") {
    if (is.null(metadata)) {
      samples <- unique(c(IBD[, "sample1"], IBD[, "sample2"]))
      vertices <- data.frame(samples)
      
    } else {
      metadata <-
        .reorder.metadata(IBD, metadata, sample.col = sample.col)
      vertices <- data.frame(metadata[, sample.col])
    }
    
    names(vertices) <- sample.col
    
    vertices
  }


.plot.IBD <-
  function(IBD,
           IBD.cutoff,
           metadata,
           label.col,
           label.palette = NULL,
           coords.cache = NULL,
           unlabelled = TRUE) {
    edgelist <- .create.edgelist(IBD)
    metadata <- .reorder.metadata(IBD, metadata)
    vertices <- .create.vertices(IBD, metadata = metadata)
    
    
    d <- edgelist[edgelist[, "fract_sites_IBD"] >= IBD.cutoff, ]
    IBD.graph <-
      graph_from_data_frame(d, directed = FALSE, vertices = vertices)
    
    
    # store metadata about the IBD graph
    if (!is.null(coords.cache)) {
      coords <- coords.cache[[as.character(IBD.cutoff)]]
    } else {
      coords <- layout_(IBD.graph, nicely())
    }
    
    if (is.null(coords)) {
      coords <- layout_(IBD.graph, nicely())
      coords.cache[[as.character(IBD.cutoff)]] <- coords
    }
    
    
    # add label and colour to the nodes
    if (is.null(label.palette)) {
      labels <- .get.labels(metadata, label.col)
      label.palette <- .generate.label.palette(labels)
      
    } else {
      labels <-
        .get.labels(metadata, label.col, label.palette = label.palette)
    }
    
    label.colours <-
      .generate.label.colours(labels, label.palette)
    
    IBD.graph <-
      set_vertex_attr(IBD.graph, "color", value = label.colours)
    
    # plot the IBD graph
    percent.cutoff <- round(IBD.cutoff * 100, digits = 2)
    
    if (unlabelled) {
      plot(
        IBD.graph,
        layout = coords,
        vertex.size = 6,
        vertex.label = NA,
        main = paste0("IBD >=", percent.cutoff, "%")
      )
    } else {
      plot(
        IBD.graph,
        layout = coords,
        vertex.size = 6,
        vertex.label.cex = 0.3,
        main = paste0("IBD >=", percent.cutoff, "%")
      )
    }
    
    
    coords.cache
  }

plot.IBD <-
  function(file,
           IBD,
           IBD.cutoffs,
           metadata,
           label.col,
           label.palette = NULL,
           onefile = TRUE,
           coords.cache = NULL,
           unlabelled = TRUE) {
    # extract relevant data to create a graph of IBDs
    edgelist <- .create.edgelist(IBD)
    metadata <- .reorder.metadata(IBD, metadata)
    
    # add label and colour to the nodes
    if (is.null(label.palette)) {
      labels <- .get.labels(metadata, label.col)
      label.palette <- .generate.label.palette(labels)
      
    } else {
      labels <-
        .get.labels(metadata, label.col, label.palette = label.palette)
    }
    
    
    if (onefile) {
      sqrt.n.plots <- ceiling(sqrt(length(IBD.cutoffs)))
      bottom.left.index <- (sqrt.n.plots - 1) * sqrt.n.plots + 1
      size <- 6 * sqrt.n.plots
      
      pdf(
        file = file,
        width = size,
        height = size
      )
      
      par(mar = c(1, 11, 1, 1) + 0.1,
          mfrow = c(sqrt.n.plots, sqrt.n.plots),
          cex.main = 1.8)
    }
    
    for (i in seq_along(IBD.cutoffs)) {
      IBD.cutoff <- IBD.cutoffs[i]
      
      if (!onefile) {
        newfile <- paste0(file, "_IBD", IBD.cutoff, ".pdf")
        pdf(file = newfile)
        par(mar = c(1, 6, 1, 1) + 0.1)
      }
      
      coords.cache <- .plot.IBD(
        IBD,
        IBD.cutoff,
        metadata,
        label.col,
        label.palette = label.palette,
        coords.cache = coords.cache,
        unlabelled = unlabelled
      )
      
      if (!onefile) {
        # add the legend about the IBD graph
        legend(
          "bottomleft",
          legend = levels(labels),
          fill = label.palette,
          cex = 0.8,
          title = label.col,
          inset = c(-0.2, 0),
          xpd = NA
        )
        dev.off()
      } else {
        if (i == bottom.left.index) {
          # add the legend about the IBD graph
          legend(
            "bottomleft",
            legend = levels(labels),
            fill = label.palette,
            cex = 1.5,
            title = label.col,
            inset = c(-0.3, 0),
            xpd = NA
          )
        }
      }
    }
    
    if (onefile)
      dev.off()
    
    invisible(coords.cache)
  }


# include -5% IBD
IBD.cutoffs <- c(0.00390625, 0.0078125, 0.015625, 0.03125, 0.0625, 0.125, 0.25, 0.5, 1) * 0.95

metadata <-
  read.csv("Popgen_db_samples_012026_filt_FINAL_fws.csv", fileEncoding = "UTF-8-BOM")



base <- metadata
for (i in 1:4) {
  copy <- base
  copy[, "wgs_id"] <- paste0(copy[, "wgs_id"], "-", i)
  metadata <- rbind(metadata, copy)
}

metadata[, "Sample"] <- metadata[, "wgs_id"]
metadata <- metadata[metadata$Year == "2024" & metadata$fws > 0.95,]



coords.cache <- vector(mode = "list", length = length(IBD.cutoffs))
names(coords.cache) <- as.character(IBD.cutoffs)

IBD.file <- "hmmIBD_new_maf0.01_out.hmm_fract.txt"
IBD <- read.delim(IBD.file)


metadata[, "district"] <- metadata[, "District"]

label.col <- "District"
IBD.plot.file <- "IBD_district.pdf"

district.palette <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#CC79A7",
  "#0072B2", "#D55E00", "#999999", "#000000", "#BBBBBB"
)
district.names <- levels(factor(metadata[, label.col]))
district.palette <- district.palette[1:length(district.names)]
names(district.palette) <- district.names

coords.cache <- plot.IBD(
  IBD.plot.file,
  IBD,
  IBD.cutoffs,
  metadata,
  label.col,
  label.palette = district.palette,
  coords.cache = coords.cache
)



######## R622I ########

metadata <-
  read.csv("Popgen_db_samples_012026_filt_FINAL_fws.csv", fileEncoding = "UTF-8-BOM")
metadata <- metadata[metadata$Year == "2024" & metadata$fws > 0.95,]

metadata_mutations <-
  read.csv("R622I_longformat.csv")

# Convert from VCF long format to wgs_id + Pfk13_R622I
metadata_mutations$Pfk13_R622I <- ifelse(
  metadata_mutations$Genotype == "0/0", "No",
  ifelse(metadata_mutations$Genotype == "1/1", "Yes",
  ifelse(metadata_mutations$Genotype == "0/1", "Het", NA))
)
colnames(metadata_mutations)[colnames(metadata_mutations) == "Sample"] <- "wgs_id"

# Add the Pfk13_R622I column to metadata
metadata <- merge(
  metadata,
  metadata_mutations[, c("wgs_id", "Pfk13_R622I")],
  by = "wgs_id",
  all.x = TRUE
)


base <- metadata
for (i in 1:4) {
  copy <- base
  copy[, "wgs_id"] <- paste0(copy[, "wgs_id"], "-", i)
  metadata <- rbind(metadata, copy)
}

metadata[, "Sample"] <- metadata[, "wgs_id"]

metadata[, "K13_R622I"] <- metadata[,"Pfk13_R622I"]

label.col <- "Pfk13_R622I"
IBD.plot.file <- "IBD_Pfk13_R622I.pdf"

coords.cache <- plot.IBD(
  IBD.plot.file,
  IBD,
  IBD.cutoffs,
  metadata,
  label.col,
  coords.cache = coords.cache
)




######## R622I Gondar only ########

metadata <-
  read.csv("Popgen_db_samples_012026_filt_FINAL_fws.csv", fileEncoding = "UTF-8-BOM")
metadata <- metadata[metadata$Year == "2024" & metadata$fws > 0.95 & metadata$District == "Gondar",]

metadata_mutations <-
  read.csv("R622I_longformat.csv")

# Convert from VCF long format to wgs_id + Pfk13_R622I
metadata_mutations$Pfk13_R622I <- ifelse(
  metadata_mutations$Genotype == "0/0", "No",
  ifelse(metadata_mutations$Genotype == "1/1", "Yes",
  ifelse(metadata_mutations$Genotype == "0/1", "Het", NA))
)
colnames(metadata_mutations)[colnames(metadata_mutations) == "Sample"] <- "wgs_id"

# Add the Pfk13_R622I column to metadata
metadata <- merge(
  metadata,
  metadata_mutations[, c("wgs_id", "Pfk13_R622I")],
  by = "wgs_id",
  all.x = TRUE
)


base <- metadata
for (i in 1:4) {
  copy <- base
  copy[, "wgs_id"] <- paste0(copy[, "wgs_id"], "-", i)
  metadata <- rbind(metadata, copy)
}

metadata[, "Sample"] <- metadata[, "wgs_id"]

label.col <- "Pfk13_R622I"
IBD.plot.file <- "IBD_Pfk13_R622I_Gondar.pdf"

IBD.file.gondar <- "hmmIBD_Gondar_maf0.01_out.hmm_fract.txt"
IBD.gondar <- read.delim(IBD.file.gondar)

coords.cache.gondar <- vector(mode = "list", length = length(IBD.cutoffs))
names(coords.cache.gondar) <- as.character(IBD.cutoffs)

coords.cache.gondar <- plot.IBD(
  IBD.plot.file,
  IBD.gondar,
  IBD.cutoffs,
  metadata,
  label.col,
  coords.cache = coords.cache.gondar
)


######## A675V ########

metadata <-
  read.csv("Popgen_db_samples_012026_filt_FINAL_fws.csv", fileEncoding = "UTF-8-BOM")
metadata <- metadata[metadata$Year == "2024" & metadata$fws > 0.95,]

metadata_mutations <-
  read.csv("A675V_longformat.csv")

# Convert from VCF long format to wgs_id + Pfk13_A675V
metadata_mutations$Pfk13_A675V <- ifelse(
  metadata_mutations$Genotype == "0/0", "No",
  ifelse(metadata_mutations$Genotype == "1/1", "Yes",
  ifelse(metadata_mutations$Genotype == "0/1", "Het", NA))
)
colnames(metadata_mutations)[colnames(metadata_mutations) == "Sample"] <- "wgs_id"

# Add the Pfk13_A675V column to metadata
metadata <- merge(
  metadata,
  metadata_mutations[, c("wgs_id", "Pfk13_A675V")],
  by = "wgs_id",
  all.x = TRUE
)



base <- metadata
for (i in 1:4) {
  copy <- base
  copy[, "wgs_id"] <- paste0(copy[, "wgs_id"], "-", i)
  metadata <- rbind(metadata, copy)
}

metadata[, "Sample"] <- metadata[, "wgs_id"]

metadata[, "K13_A675V"] <- metadata[,"Pfk13_A675V"]

label.col <- "Pfk13_A675V"
IBD.plot.file <- "IBD_Pfk13_A675V.pdf"

coords.cache <- plot.IBD(
  IBD.plot.file,
  IBD,
  IBD.cutoffs,
  metadata,
  label.col,
  coords.cache = coords.cache
)


######## HRP deletion ########

metadata <-
  read.csv("Popgen_db_samples_012026_filt_FINAL_fws.csv", fileEncoding = "UTF-8-BOM")
metadata <- metadata[metadata$Year == "2024" & metadata$fws > 0.95,]

metadata_hrp <-
  read.csv("sample_meta_with_hrp_status_updated.csv")

colnames(metadata_hrp)[colnames(metadata_hrp) == "Sample"] <- "wgs_id"

# Add the hrp_status column to metadata
metadata <- merge(
  metadata,
  metadata_hrp[, c("wgs_id", "hrp_status")],
  by = "wgs_id",
  all.x = TRUE
)



base <- metadata
for (i in 1:4) {
  copy <- base
  copy[, "wgs_id"] <- paste0(copy[, "wgs_id"], "-", i)
  metadata <- rbind(metadata, copy)
}

metadata[, "Sample"] <- metadata[, "wgs_id"]

# Recode hrp_status to proper nomenclature
metadata$hrp_label <- ifelse(
  metadata$hrp_status == "wild", "Pfhrp2+/Pfhrp3+",
  ifelse(metadata$hrp_status == "hrp2_only", "Pfhrp2-/Pfhrp3+",
  ifelse(metadata$hrp_status == "hrp3_only", "Pfhrp2+/Pfhrp3-",
  ifelse(metadata$hrp_status == "mutant", "Pfhrp2-/Pfhrp3-", NA)))
)

label.col <- "hrp_label"
IBD.plot.file <- "IBD_hrp_status.pdf"

coords.cache <- plot.IBD(
  IBD.plot.file,
  IBD,
  IBD.cutoffs,
  metadata,
  label.col,
  coords.cache = coords.cache
)



######## District + HRP deletion (colour + shape) ########

metadata <-
  read.csv("Popgen_db_samples_012026_filt_FINAL_fws.csv", fileEncoding = "UTF-8-BOM")
metadata <- metadata[metadata$Year == "2024" & metadata$fws > 0.95,]

metadata_hrp <-
  read.csv("sample_meta_with_hrp_status_updated.csv")
colnames(metadata_hrp)[colnames(metadata_hrp) == "Sample"] <- "wgs_id"

metadata <- merge(
  metadata,
  metadata_hrp[, c("wgs_id", "hrp_status")],
  by = "wgs_id",
  all.x = TRUE
)


base <- metadata
for (i in 1:4) {
  copy <- base
  copy[, "wgs_id"] <- paste0(copy[, "wgs_id"], "-", i)
  metadata <- rbind(metadata, copy)
}
metadata[, "Sample"] <- metadata[, "wgs_id"]

# Map HRP status to shapes
hrp.shape.map <- c(
  "wild"      = "circle",
  "hrp2_only" = "square",
  "hrp3_only" = "triangle",
  "mutant"    = "diamond"
)
metadata$hrp_shape <- hrp.shape.map[metadata$hrp_status]
metadata$hrp_shape[is.na(metadata$hrp_shape)] <- "circle"

IBD.plot.file <- "IBD_district_hrp.pdf"

# Reuse district palette
district.labels <- .get.labels(metadata, "District")
if (!exists("district.palette")) {
  dp <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#CC79A7",
          "#0072B2", "#D55E00", "#999999", "#000000", "#BBBBBB")
  district.names <- levels(district.labels)
  district.palette <- dp[1:length(district.names)]
  names(district.palette) <- district.names
}
district.colours <- district.palette[district.labels]

# Create the plot
edgelist <- .create.edgelist(IBD)
metadata <- .reorder.metadata(IBD, metadata)
vertices <- .create.vertices(IBD, metadata = metadata)

# Recalculate colours and shapes after reorder
district.labels <- .get.labels(metadata, "District")
district.colours <- district.palette[district.labels]
vertex.shapes <- metadata$hrp_shape

sqrt.n.plots <- ceiling(sqrt(length(IBD.cutoffs)))
bottom.left.index <- (sqrt.n.plots - 1) * sqrt.n.plots + 1
size <- 6 * sqrt.n.plots

pdf(file = IBD.plot.file, width = size, height = size)
par(mar = c(1, 11, 1, 1) + 0.1,
    mfrow = c(sqrt.n.plots, sqrt.n.plots),
    cex.main = 1.8)

for (i in seq_along(IBD.cutoffs)) {
  IBD.cutoff <- IBD.cutoffs[i]
  d <- edgelist[edgelist[, "fract_sites_IBD"] >= IBD.cutoff, ]
  IBD.graph <- graph_from_data_frame(d, directed = FALSE, vertices = vertices)

  coords <- coords.cache[[as.character(IBD.cutoff)]]
  if (is.null(coords)) {
    coords <- layout_(IBD.graph, nicely())
    coords.cache[[as.character(IBD.cutoff)]] <- coords
  }

  percent.cutoff <- round(IBD.cutoff * 100, digits = 2)

  # Map HRP status to pch values
  pch.map <- c("circle" = 21, "square" = 22, "triangle" = 24, "diamond" = 23)
  vertex.pch <- pch.map[vertex.shapes]
  vertex.pch[is.na(vertex.pch)] <- 21

  # Plot graph with invisible vertices
  plot(
    IBD.graph,
    layout = coords,
    vertex.size = 6,
    vertex.shape = "none",
    vertex.label = NA,
    main = paste0("IBD >=", percent.cutoff, "%"),
    rescale = TRUE
  )

  # Get rescaled coordinates matching what igraph plotted
  layout.norm <- norm_coords(coords, xmin = -1, xmax = 1, ymin = -1, ymax = 1)

  # Draw vertices manually with correct shapes and colours
  points(layout.norm[, 1], layout.norm[, 2],
         pch = vertex.pch, bg = district.colours, col = "black",
         cex = 2, lwd = 0.5)

  if (i == bottom.left.index) {
    # District colour legend
    legend(
      "bottomleft",
      legend = names(district.palette),
      fill = district.palette,
      cex = 1.2,
      title = "District",
      inset = c(-0.3, 0),
      xpd = NA
    )
    # HRP shape legend
    hrp.labels <- c(
      "wild"      = expression(italic("Pfhrp2") * "+/" * italic("Pfhrp3") * "+"),
      "hrp2_only" = expression(italic("Pfhrp2") * "-/" * italic("Pfhrp3") * "+"),
      "hrp3_only" = expression(italic("Pfhrp2") * "+/" * italic("Pfhrp3") * "-"),
      "mutant"    = expression(italic("Pfhrp2") * "-/" * italic("Pfhrp3") * "-")
    )
    legend(
      "bottomright",
      legend = hrp.labels,
      pch = c(21, 22, 24, 23),  # circle, square, triangle, diamond
      pt.bg = "grey50",
      cex = 1.2,
      title = "HRP status",
      inset = c(-0.1, 0),
      xpd = NA
    )
  }
}

dev.off()


######## Boxplots: pairwise IBD among mutants vs wildtypes per district ########

library(ggplot2)
library(dplyr)

# Function to build a dataframe of pairwise IBD annotated with mutation status and district
.build.ibd.by.status <- function(IBD, metadata, mutation.col) {
  sample.district <- setNames(metadata$District, metadata$Sample)
  sample.status <- setNames(metadata[, mutation.col], metadata$Sample)

  s1.status <- sample.status[IBD$sample1]
  s2.status <- sample.status[IBD$sample2]
  s1.district <- sample.district[IBD$sample1]
  s2.district <- sample.district[IBD$sample2]

  same.district <- !is.na(s1.district) & !is.na(s2.district) & s1.district == s2.district
  both.known <- !is.na(s1.status) & !is.na(s2.status)
  keep <- same.district & both.known

  pair.status <- ifelse(
    s1.status[keep] == "Yes" & s2.status[keep] == "Yes", "Mutant-Mutant",
    ifelse(
      s1.status[keep] == "No" & s2.status[keep] == "No", "WT-WT",
      "Mixed"
    )
  )

  data.frame(
    fract_sites_IBD = IBD$fract_sites_IBD[keep],
    district = s1.district[keep],
    pair_status = factor(pair.status, levels = c("WT-WT", "Mutant-Mutant"),
                         labels = c("WT", "Mutant")),
    stringsAsFactors = FALSE
  )
}

.plot.ibd.boxplot <- function(ibd.df, title, filename, metadata, mutation.col) {
  
  library(dplyr)
  library(ggplot2)
  
  # Keep only WT vs Mutant pairs
  ibd.df <- ibd.df[ibd.df$pair_status %in% c("WT", "Mutant"),]

  # ---------------------------
  # Sample counts (per district)
  # ---------------------------
  meta.sub <- metadata[!is.na(metadata[, mutation.col]) & !is.na(metadata$District),]
  meta.sub$pair_status <- ifelse(meta.sub[, mutation.col] == "Yes", "Mutant",
                           ifelse(meta.sub[, mutation.col] == "No", "WT", NA))
  meta.sub <- meta.sub[!is.na(meta.sub$pair_status),]

  counts <- aggregate(Sample ~ pair_status + District, data = meta.sub, FUN = length)
  colnames(counts) <- c("pair_status", "district", "n")
  counts$label <- paste0("n=", counts$n)

  # ---------------------------
  # Y positions
  # ---------------------------
  ymax <- aggregate(fract_sites_IBD ~ district, data = ibd.df, FUN = max)
  colnames(ymax)[2] <- "ymax"

  counts <- merge(counts, ymax, by = "district")
  counts$y_n <- counts$ymax * 1.05   # n labels

  # ---------------------------
  # P-values per district
  # ---------------------------
  pvals <- do.call(rbind, lapply(split(ibd.df, ibd.df$district), function(df) {
    if (length(unique(df$pair_status)) < 2) {
      return(data.frame(district = unique(df$district), p_value = NA))
    }
    data.frame(
      district = unique(df$district),
      p_value = wilcox.test(fract_sites_IBD ~ pair_status, data = df)$p.value
    )
  }))

  pvals <- merge(pvals, ymax, by = "district")
  pvals <- pvals[!is.na(pvals$p_value),]

  # Format p-values as significance stars
  pvals$label <- ifelse(pvals$p_value < 0.001, "***",
                 ifelse(pvals$p_value < 0.01, "**",
                 ifelse(pvals$p_value < 0.05, "*", "ns")))
  pvals$y_line <- pvals$ymax * 1.3
  pvals$y_p <- pvals$ymax * 1.4

  # ---------------------------
  # Plot
  # ---------------------------
  p <- ggplot(ibd.df, aes(x = pair_status, y = fract_sites_IBD,
                          fill = pair_status, colour = pair_status)) +
    
    geom_violin(alpha = 0.3, scale = "width") +
    geom_jitter(size = 1, alpha = 0.5, width = 0.15) +
    geom_boxplot(width = 0.15, fill = "white", alpha = 0.7,
                 outlier.shape = NA, colour = "black") +
    
    # n labels
    geom_text(
      data = counts,
      aes(x = pair_status, y = y_n, label = label),
      vjust = 0,
      size = 3,
      colour = "black"
    ) +
    
    # horizontal comparison line
    geom_segment(
      data = pvals,
      aes(x = 1, xend = 2, y = y_line, yend = y_line),
      inherit.aes = FALSE,
      linewidth = 0.4
    ) +
    
    # p-values
    geom_text(
      data = pvals,
      aes(x = 1.5, y = y_p, label = label),
      inherit.aes = FALSE,
      size = 3
    ) +
    
    facet_wrap(~district, scales = "free_y") +
    
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.30))) +
    
    scale_fill_manual(values = c("WT" = "#A6CEE3", "Mutant" = "#E31A1C")) +
    scale_colour_manual(values = c("WT" = "#6BAED6", "Mutant" = "#CB181D")) +
    
    labs(title = title,
         y = "Pairwise IBD fraction",
         x = "",
         fill = "Pair status",
         colour = "Pair status") +
    
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.background = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA)
    )

  ggsave(filename, p, width = 10, height = 6)
}




# --- R622I boxplot ---
metadata_r622i <-
  read.csv("Popgen_db_samples_012026_filt_FINAL_fws.csv", fileEncoding = "UTF-8-BOM")
metadata_r622i <- metadata_r622i[metadata_r622i$Year == "2024" & metadata_r622i$fws > 0.95,]
mut <- read.csv("R622I_longformat.csv")
mut$Pfk13_R622I <- ifelse(
  mut$Genotype == "0/0", "No",
  ifelse(mut$Genotype == "1/1", "Yes",
  ifelse(mut$Genotype == "0/1", "Het", NA))
)
colnames(mut)[colnames(mut) == "Sample"] <- "wgs_id"
metadata_r622i <- merge(metadata_r622i, mut[, c("wgs_id", "Pfk13_R622I")], by = "wgs_id", all.x = TRUE)
metadata_r622i$Sample <- metadata_r622i$wgs_id

ibd.df <- .build.ibd.by.status(IBD, metadata_r622i, "Pfk13_R622I")
.plot.ibd.boxplot(ibd.df, "Pairwise IBD: Pfk13 R622I mutant vs WT by district", "IBD_boxplot_R622I.pdf", metadata_r622i, "Pfk13_R622I")


# --- R622I Gondar boxplot ---
metadata_r622i_gondar <-
  read.csv("Popgen_db_samples_012026_filt_FINAL_fws.csv", fileEncoding = "UTF-8-BOM")
metadata_r622i_gondar <- metadata_r622i_gondar[metadata_r622i_gondar$Year == "2024" & metadata_r622i_gondar$fws > 0.95 & metadata_r622i_gondar$District == "Gondar",]
mut <- read.csv("R622I_longformat.csv")
mut$Pfk13_R622I <- ifelse(
  mut$Genotype == "0/0", "No",
  ifelse(mut$Genotype == "1/1", "Yes",
  ifelse(mut$Genotype == "0/1", "Het", NA))
)
colnames(mut)[colnames(mut) == "Sample"] <- "wgs_id"
metadata_r622i_gondar <- merge(metadata_r622i_gondar, mut[, c("wgs_id", "Pfk13_R622I")], by = "wgs_id", all.x = TRUE)
metadata_r622i_gondar$Sample <- metadata_r622i_gondar$wgs_id

ibd.df <- .build.ibd.by.status(IBD.gondar, metadata_r622i_gondar, "Pfk13_R622I")
.plot.ibd.boxplot(ibd.df, "Pairwise IBD: Pfk13 R622I mutant vs WT (Gondar)", "IBD_boxplot_R622I_Gondar.pdf", metadata_r622i_gondar, "Pfk13_R622I")


# --- A675V boxplot ---
metadata_a675v <-
  read.csv("Popgen_db_samples_012026_filt_FINAL_fws.csv", fileEncoding = "UTF-8-BOM")
metadata_a675v <- metadata_a675v[metadata_a675v$Year == "2024" & metadata_a675v$fws > 0.95,]
mut <- read.csv("A675V_longformat.csv")
mut$Pfk13_A675V <- ifelse(
  mut$Genotype == "0/0", "No",
  ifelse(mut$Genotype == "1/1", "Yes",
  ifelse(mut$Genotype == "0/1", "Het", NA))
)
colnames(mut)[colnames(mut) == "Sample"] <- "wgs_id"
metadata_a675v <- merge(metadata_a675v, mut[, c("wgs_id", "Pfk13_A675V")], by = "wgs_id", all.x = TRUE)
metadata_a675v$Sample <- metadata_a675v$wgs_id

ibd.df <- .build.ibd.by.status(IBD, metadata_a675v, "Pfk13_A675V")
.plot.ibd.boxplot(ibd.df, "Pairwise IBD: Pfk13 A675V mutant vs WT by district", "IBD_boxplot_A675V.pdf", metadata_a675v, "Pfk13_A675V")


# --- HRP deletion boxplot ---
metadata_hrp_bp <-
  read.csv("Popgen_db_samples_012026_filt_FINAL_fws.csv", fileEncoding = "UTF-8-BOM")
metadata_hrp_bp <- metadata_hrp_bp[metadata_hrp_bp$Year == "2024" & metadata_hrp_bp$fws > 0.95,]
hrp <- read.csv("sample_meta_with_hrp_status_updated.csv")
colnames(hrp)[colnames(hrp) == "Sample"] <- "wgs_id"
hrp$hrp_mutant <- ifelse(hrp$hrp_status == "wild", "No", "Yes")
metadata_hrp_bp <- merge(metadata_hrp_bp, hrp[, c("wgs_id", "hrp_mutant")], by = "wgs_id", all.x = TRUE)
metadata_hrp_bp$Sample <- metadata_hrp_bp$wgs_id

ibd.df <- .build.ibd.by.status(IBD, metadata_hrp_bp, "hrp_mutant")
.plot.ibd.boxplot(ibd.df, "Pairwise IBD: HRP deletion vs WT by district", "IBD_boxplot_hrp.pdf", metadata_hrp_bp, "hrp_mutant")


######## Gondar + Eritrean samples: District colour + R622I shape  ########

IBD.region.r622i <- read.delim("hmmIBD_region_R622I_maf0.01_out.hmm_fract.txt")

metadata_region <- read.csv("Popgen_db_samples_012026_filt_FINAL_fws.csv", fileEncoding = "UTF-8-BOM")
metadata_region <- metadata_region[
  metadata_region$fws > 0.95 &
  (metadata_region$Country == "Eritrea" | metadata_region$District == "Gondar"), ]

mut_r622i <- read.csv("R622I_longformat.csv")

# Normalize pipe vs slash genotype notation
mut_r622i$Genotype <- with(mut_r622i, ifelse(
  Genotype == "1|1", "1/1", ifelse(
  Genotype == "0|1" | Genotype == "1|0", "0/1", ifelse(
  Genotype == "0|0", "0/0", Genotype))))

# Deduplicate: one row per sample
genotype_rank_vec <- c("1/1" = 3, "0/1" = 2, "0/0" = 1)
mut_r622i$rank <- genotype_rank_vec[mut_r622i$Genotype]
mut_r622i$rank[is.na(mut_r622i$rank)] <- 0
mut_r622i <- mut_r622i[order(mut_r622i$Sample, -mut_r622i$rank), ]
mut_r622i <- mut_r622i[!duplicated(mut_r622i$Sample), ]

mut_r622i$Pfk13_R622I <- ifelse(
  mut_r622i$Genotype == "1/1", "Yes",
  ifelse(mut_r622i$Genotype == "0/1", "Het",
  ifelse(mut_r622i$Genotype == "0/0", "No", NA)))
colnames(mut_r622i)[colnames(mut_r622i) == "Sample"] <- "wgs_id"

metadata_region <- merge(
  metadata_region,
  mut_r622i[, c("wgs_id", "Pfk13_R622I")],
  by = "wgs_id", all.x = TRUE
)

base <- metadata_region
for (i in 1:4) {
  copy <- base
  copy[, "wgs_id"] <- paste0(copy[, "wgs_id"], "-", i)
  metadata_region <- rbind(metadata_region, copy)
}
metadata_region[, "Sample"] <- metadata_region[, "wgs_id"]

# Filter IBD to only pairs where both samples are in filtered metadata
keep.samples <- metadata_region$Sample
IBD.region.r622i <- IBD.region.r622i[
  IBD.region.r622i$sample1 %in% keep.samples &
  IBD.region.r622i$sample2 %in% keep.samples, ]

# Re-assert District from Admin.level.1 for any Eritrean samples that lost it after merge/rbind
metadata_region$District <- ifelse(
  metadata_region$Country == "Eritrea" & (is.na(metadata_region$District) | metadata_region$District == ""),
  metadata_region$Admin.level.1,
  metadata_region$District
)

# Drop any remaining rows with NA District (safety net)
metadata_region <- metadata_region[!is.na(metadata_region$District) & metadata_region$District != "", ]

# Map R622I to shapes
r622i.shape.map <- c("No" = "circle", "Yes" = "triangle", "Het" = "diamond")
metadata_region$r622i_shape <- r622i.shape.map[metadata_region$Pfk13_R622I]
metadata_region$r622i_shape[is.na(metadata_region$r622i_shape)] <- "circle"

# Also re-filter IBD after dropping any NA-district rows
keep.samples <- metadata_region$Sample
IBD.region.r622i <- IBD.region.r622i[
  IBD.region.r622i$sample1 %in% keep.samples &
  IBD.region.r622i$sample2 %in% keep.samples, ]

# District palette for this dataset
district.labels.region <- .get.labels(metadata_region, "District")
district.names.region  <- levels(district.labels.region)
dp <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#CC79A7",
        "#0072B2", "#D55E00", "#999999", "#000000", "#BBBBBB",
        "#44AA99", "#882255", "#DDCC77", "#117733", "#332288")
district.palette.region <- dp[seq_along(district.names.region)]
names(district.palette.region) <- district.names.region

# Build graph objects
edgelist.region  <- .create.edgelist(IBD.region.r622i)
metadata_region  <- .reorder.metadata(IBD.region.r622i, metadata_region)
vertices.region  <- .create.vertices(IBD.region.r622i, metadata = metadata_region)

district.labels.region  <- .get.labels(metadata_region, "District")
district.colours.region <- district.palette.region[district.labels.region]
vertex.shapes.region    <- metadata_region$r622i_shape

pch.map <- c("circle" = 21, "square" = 22, "triangle" = 24, "diamond" = 23)

coords.cache.region <- vector(mode = "list", length = length(IBD.cutoffs))
names(coords.cache.region) <- as.character(IBD.cutoffs)

sqrt.n.plots      <- ceiling(sqrt(length(IBD.cutoffs)))
bottom.left.index <- (sqrt.n.plots - 1) * sqrt.n.plots + 1
size              <- 6 * sqrt.n.plots

pdf(file = "IBD_district_R622I_region.pdf", width = size, height = size)
par(mar = c(1, 11, 1, 1) + 0.1,
    mfrow = c(sqrt.n.plots, sqrt.n.plots),
    cex.main = 1.8)

for (i in seq_along(IBD.cutoffs)) {
  IBD.cutoff <- IBD.cutoffs[i]
  d         <- edgelist.region[edgelist.region[, "fract_sites_IBD"] >= IBD.cutoff, ]
  IBD.graph <- graph_from_data_frame(d, directed = FALSE, vertices = vertices.region)

  coords <- coords.cache.region[[as.character(IBD.cutoff)]]
  if (is.null(coords)) {
    coords <- layout_(IBD.graph, nicely())
    coords.cache.region[[as.character(IBD.cutoff)]] <- coords
  }

  percent.cutoff <- round(IBD.cutoff * 100, digits = 2)
  vertex.pch <- pch.map[vertex.shapes.region]
  vertex.pch[is.na(vertex.pch)] <- 21

  plot(
    IBD.graph,
    layout      = coords,
    vertex.size  = 6,
    vertex.shape = "none",
    vertex.label = NA,
    main         = paste0("IBD >=", percent.cutoff, "%"),
    rescale      = TRUE
  )

  layout.norm <- norm_coords(coords, xmin = -1, xmax = 1, ymin = -1, ymax = 1)
  points(layout.norm[, 1], layout.norm[, 2],
         pch = vertex.pch, bg = district.colours.region, col = "black",
         cex = 2, lwd = 0.5)

  if (i == bottom.left.index) {
    legend(
      "bottomleft",
      legend = names(district.palette.region),
      fill   = district.palette.region,
      cex    = 1.2,
      title  = "District",
      inset  = c(-0.3, 0),
      xpd    = NA
    )
    legend(
      "bottomright",
      legend = c("WT", "Heterozygous", "Mutant"),
      pch    = c(21, 23, 24),
      pt.bg  = "grey50",
      cex    = 1.2,
      title  = "PfK13 R622I",
      inset  = c(-0.1, 0),
      xpd    = NA
    )
  }
}

dev.off()


# --- Regional R622I boxplot per district ---
metadata_r622i_region_bp <- read.csv("Popgen_db_samples_012026_filt_FINAL_fws.csv", fileEncoding = "UTF-8-BOM")
metadata_r622i_region_bp <- metadata_r622i_region_bp[metadata_r622i_region_bp$fws > 0.95, ]
metadata_r622i_region_bp <- merge(
  metadata_r622i_region_bp,
  mut_r622i[, c("wgs_id", "Pfk13_R622I")],
  by = "wgs_id", all.x = TRUE
)
metadata_r622i_region_bp$Sample <- metadata_r622i_region_bp$wgs_id

ibd.df.region <- .build.ibd.by.status(IBD.region.r622i, metadata_r622i_region_bp, "Pfk13_R622I")
.plot.ibd.boxplot(
  ibd.df.region,
  "Pairwise IBD: PfK13 R622I mutant vs WT by district (regional)",
  "IBD_boxplot_R622I_region.pdf",
  metadata_r622i_region_bp,
  "Pfk13_R622I"
)


######## Gondar only: R622I + hrp2/3 deletion co-carriage ########

metadata_gondar <-
  read.csv("Popgen_db_samples_012026_filt_FINAL_fws.csv", fileEncoding = "UTF-8-BOM")
metadata_gondar <- metadata_gondar[
  metadata_gondar$Year == "2024" &
  metadata_gondar$fws > 0.95 &
  metadata_gondar$District == "Gondar", ]

# --- R622I genotype and hrp2/3 deletion status, from the combined table ---
# (one row per sample, so no deduplication is needed)
combined_gondar <- read.csv("sample_meta_R622I_combined.csv")
combined_gondar$Pfk13_R622I <- ifelse(
  combined_gondar$Genotype %in% c("0/0", "0|0"), "No",
  ifelse(combined_gondar$Genotype %in% c("1/1", "1|1"), "Yes",
  ifelse(combined_gondar$Genotype %in% c("0/1", "0|1", "1|0"), "Het", NA))
)
colnames(combined_gondar)[colnames(combined_gondar) == "Sample"] <- "wgs_id"

metadata_gondar <- merge(
  metadata_gondar,
  combined_gondar[, c("wgs_id", "Pfk13_R622I", "hrp_status")],
  by = "wgs_id", all.x = TRUE
)

# Carrier definitions. Het R622I calls are counted as carriers (none are present
# among the Gondar samples in the IBD matrix, but keep the rule explicit); any
# hrp_status other than "wild" is a hrp2 and/or hrp3 deletion.
r622i.carrier <- metadata_gondar$Pfk13_R622I %in% c("Yes", "Het")
hrp.carrier   <- !is.na(metadata_gondar$hrp_status) & metadata_gondar$hrp_status != "wild"
genotyped     <- !is.na(metadata_gondar$Pfk13_R622I) & !is.na(metadata_gondar$hrp_status)

metadata_gondar$combo_status <- ifelse(
  !genotyped, "Not genotyped",
  ifelse(r622i.carrier & hrp.carrier, "R622I + hrp2/3 deletion",
  ifelse(r622i.carrier, "R622I only",
  ifelse(hrp.carrier, "hrp2/3 deletion only", "Neither"))))

combo.levels <- c("Neither", "R622I only", "hrp2/3 deletion only",
                  "R622I + hrp2/3 deletion", "Not genotyped")
combo.palette <- c(
  "Neither"                 = "#BBBBBB",
  "R622I only"              = "#E69F00",
  "hrp2/3 deletion only"    = "#56B4E9",
  "R622I + hrp2/3 deletion" = "#D7191C",
  "Not genotyped"           = "#FFFFFF"
)

# Shape encoding for the second plot: shape = hrp status, colour = R622I status
metadata_gondar$hrp_shape <- hrp.shape.map[metadata_gondar$hrp_status]
metadata_gondar$hrp_shape[is.na(metadata_gondar$hrp_shape)] <- "circle"

r622i.palette <- c("No" = "#A6CEE3", "Het" = "#FDBF6F", "Yes" = "#D7191C")

# Duplicated ids (-1 .. -4) mirror the convention used elsewhere in this script
base <- metadata_gondar
for (i in 1:4) {
  copy <- base
  copy[, "wgs_id"] <- paste0(copy[, "wgs_id"], "-", i)
  metadata_gondar <- rbind(metadata_gondar, copy)
}
metadata_gondar[, "Sample"] <- metadata_gondar[, "wgs_id"]

# Build graph objects (IBD.gondar / coords.cache.gondar come from the
# "R622I Gondar only" section above)
edgelist.gondar <- .create.edgelist(IBD.gondar)
metadata_gondar <- .reorder.metadata(IBD.gondar, metadata_gondar)
vertices.gondar <- .create.vertices(IBD.gondar, metadata = metadata_gondar)

# Recompute per-vertex aesthetics after the reorder
combo.labels  <- factor(metadata_gondar$combo_status, levels = combo.levels)
combo.colours <- combo.palette[combo.labels]
combo.pch     <- ifelse(metadata_gondar$combo_status == "R622I + hrp2/3 deletion", 23, 21)
combo.cex     <- ifelse(metadata_gondar$combo_status == "R622I + hrp2/3 deletion", 2.8, 2)
combo.lwd     <- ifelse(metadata_gondar$combo_status == "R622I + hrp2/3 deletion", 1.5, 0.5)

r622i.colours.gondar <- r622i.palette[metadata_gondar$Pfk13_R622I]
r622i.colours.gondar[is.na(r622i.colours.gondar)] <- "#FFFFFF"
vertex.pch.gondar <- pch.map[metadata_gondar$hrp_shape]
vertex.pch.gondar[is.na(vertex.pch.gondar)] <- 21

# Counts per category, for the legend
combo.counts <- table(factor(metadata_gondar$combo_status, levels = combo.levels))
combo.keep   <- combo.levels[combo.counts > 0]

sqrt.n.plots      <- ceiling(sqrt(length(IBD.cutoffs)))
bottom.left.index <- (sqrt.n.plots - 1) * sqrt.n.plots + 1
size              <- 6 * sqrt.n.plots

.gondar.coords <- function(IBD.graph, IBD.cutoff) {
  coords <- coords.cache.gondar[[as.character(IBD.cutoff)]]
  # only reuse a cached layout if it matches this graph's vertex set
  if (is.null(coords) || nrow(coords) != vcount(IBD.graph)) {
    coords <- layout_(IBD.graph, nicely())
    coords.cache.gondar[[as.character(IBD.cutoff)]] <<- coords
  }
  coords
}


# --- Plot 1: single colour per co-carriage category ---
pdf(file = "IBD_Gondar_R622I_hrp.pdf", width = size, height = size)
par(mar = c(1, 11, 1, 1) + 0.1,
    mfrow = c(sqrt.n.plots, sqrt.n.plots),
    cex.main = 1.8)

for (i in seq_along(IBD.cutoffs)) {
  IBD.cutoff <- IBD.cutoffs[i]
  d <- edgelist.gondar[edgelist.gondar[, "fract_sites_IBD"] >= IBD.cutoff, ]
  IBD.graph <- graph_from_data_frame(d, directed = FALSE, vertices = vertices.gondar)

  coords <- .gondar.coords(IBD.graph, IBD.cutoff)
  percent.cutoff <- round(IBD.cutoff * 100, digits = 2)

  plot(
    IBD.graph,
    layout       = coords,
    vertex.size  = 6,
    vertex.shape = "none",
    vertex.label = NA,
    main         = paste0("IBD >=", percent.cutoff, "%"),
    rescale      = TRUE
  )

  layout.norm <- norm_coords(coords, xmin = -1, xmax = 1, ymin = -1, ymax = 1)

  # draw co-carriers last so they sit on top
  draw.order <- order(metadata_gondar$combo_status == "R622I + hrp2/3 deletion")
  points(layout.norm[draw.order, 1], layout.norm[draw.order, 2],
         pch = combo.pch[draw.order], bg = combo.colours[draw.order],
         col = "black", cex = combo.cex[draw.order], lwd = combo.lwd[draw.order])

  if (i == bottom.left.index) {
    legend(
      "bottomleft",
      legend = paste0(combo.keep, " (n=", combo.counts[combo.keep], ")"),
      pch    = ifelse(combo.keep == "R622I + hrp2/3 deletion", 23, 21),
      pt.bg  = combo.palette[combo.keep],
      pt.cex = ifelse(combo.keep == "R622I + hrp2/3 deletion", 2, 1.5),
      cex    = 1.2,
      title  = "PfK13 R622I / hrp2-3 status",
      inset  = c(-0.25, 0),
      xpd    = NA
    )
  }
}

dev.off()


# --- Plot 2: colour = R622I status, shape = hrp2/3 status ---
pdf(file = "IBD_Gondar_R622I_hrp_shapes.pdf", width = size, height = size)
par(mar = c(1, 11, 1, 1) + 0.1,
    mfrow = c(sqrt.n.plots, sqrt.n.plots),
    cex.main = 1.8)

for (i in seq_along(IBD.cutoffs)) {
  IBD.cutoff <- IBD.cutoffs[i]
  d <- edgelist.gondar[edgelist.gondar[, "fract_sites_IBD"] >= IBD.cutoff, ]
  IBD.graph <- graph_from_data_frame(d, directed = FALSE, vertices = vertices.gondar)

  coords <- .gondar.coords(IBD.graph, IBD.cutoff)
  percent.cutoff <- round(IBD.cutoff * 100, digits = 2)

  plot(
    IBD.graph,
    layout       = coords,
    vertex.size  = 6,
    vertex.shape = "none",
    vertex.label = NA,
    main         = paste0("IBD >=", percent.cutoff, "%"),
    rescale      = TRUE
  )

  layout.norm <- norm_coords(coords, xmin = -1, xmax = 1, ymin = -1, ymax = 1)
  points(layout.norm[, 1], layout.norm[, 2],
         pch = vertex.pch.gondar, bg = r622i.colours.gondar,
         col = "black", cex = 2, lwd = 0.5)

  if (i == bottom.left.index) {
    legend(
      "bottomleft",
      legend = c("WT", "Heterozygous", "Mutant"),
      pch    = 21,
      pt.bg  = r622i.palette[c("No", "Het", "Yes")],
      pt.cex = 1.5,
      cex    = 1.2,
      title  = "PfK13 R622I",
      inset  = c(-0.3, 0),
      xpd    = NA
    )
    legend(
      "bottomright",
      legend = c(
        expression(italic("Pfhrp2") * "+/" * italic("Pfhrp3") * "+"),
        expression(italic("Pfhrp2") * "-/" * italic("Pfhrp3") * "+"),
        expression(italic("Pfhrp2") * "+/" * italic("Pfhrp3") * "-"),
        expression(italic("Pfhrp2") * "-/" * italic("Pfhrp3") * "-")
      ),
      pch    = c(21, 22, 24, 23),
      pt.bg  = "grey50",
      pt.cex = 1.5,
      cex    = 1.2,
      title  = "HRP status",
      inset  = c(-0.1, 0),
      xpd    = NA
    )
  }
}

dev.off()


# --- Co-carriage cross-tabulation (one row per sample, printed to console) ---
gondar.unique <- metadata_gondar[!grepl("-[1-4]$", metadata_gondar$Sample), ]
cat("\nGondar 2024 (fws > 0.95): PfK13 R622I vs hrp2/3 status\n")
print(table(R622I = gondar.unique$Pfk13_R622I,
            hrp   = gondar.unique$hrp_status, useNA = "ifany"))


######## Gambella + Mizan: PfPX1 PIN haplotype ########

# No dedicated hmmIBD run is needed: every Gambella/Mizan sample and all of their
# pairwise comparisons are present in the all-Ethiopia output, so the subset is
# taken from there. Note that the IBD estimates therefore use allele frequencies
# pooled across all Ethiopian districts.

px1.districts <- c("Gambella", "Mizan")

metadata_px1 <-
  read.csv("Popgen_db_samples_012026_filt_FINAL_fws.csv", fileEncoding = "UTF-8-BOM")
metadata_px1 <- metadata_px1[
  metadata_px1$Year == "2024" &
  metadata_px1$fws > 0.95 &
  metadata_px1$District %in% px1.districts, ]

# --- PfPX1 genotypes: one row per sample, one column per codon ---
# Calls are "WT", "L1222P"/"M1701I"/"D1705N", or NA (1222 L>P, 1701 M>I,
# 1705 D>N). read.delim turns the numeric headers into X1222 / X1701 / X1705.
# This file has no heterozygous calls, so there is no "Mixed call" class.
px1.wide <- read.delim("px1_genotype_result.tsv")

px1.idx <- match(metadata_px1$wgs_id, px1.wide$Sample)
metadata_px1$px1_1222 <- px1.wide$X1222[px1.idx]
metadata_px1$px1_1701 <- px1.wide$X1701[px1.idx]
metadata_px1$px1_1705 <- px1.wide$X1705[px1.idx]

px1.codons <- metadata_px1[, c("px1_1222", "px1_1701", "px1_1705")]

n.called <- rowSums(!is.na(px1.codons))
n.mutant <- rowSums(!is.na(px1.codons) & px1.codons != "WT")

metadata_px1$px1_haplotype <- ifelse(
  n.called < 3, "Incomplete call",
  ifelse(n.mutant == 3, "PIN (triple mutant)",
  ifelse(n.mutant == 0, "LMD (wild type)", "Partial (1-2 mutations)")))

# Amino-acid haplotype string (e.g. "LMD", "PIN", "LID"), for the console summary
px1.wild <- c("L", "M", "D")              # LMD haplotype; mutant is PIN
metadata_px1$px1_hap_string <- apply(
  mapply(function(call, ref)
    ifelse(is.na(call), "?", ifelse(call == "WT", ref, substr(call, nchar(call), nchar(call)))),
    px1.codons, px1.wild),
  1, paste, collapse = "")

px1.levels <- c("LMD (wild type)", "Partial (1-2 mutations)", "PIN (triple mutant)",
                "Incomplete call")
px1.palette <- c(
  "LMD (wild type)"         = "#A6CEE3",
  "Partial (1-2 mutations)" = "#E69F00",
  "PIN (triple mutant)"     = "#D7191C",
  "Incomplete call"         = "#FFFFFF"
)

# District is encoded by shape
px1.shape.map <- c("Gambella" = "circle", "Mizan" = "square")
metadata_px1$px1_shape <- px1.shape.map[metadata_px1$District]
metadata_px1$px1_shape[is.na(metadata_px1$px1_shape)] <- "circle"

base <- metadata_px1
for (i in 1:4) {
  copy <- base
  copy[, "wgs_id"] <- paste0(copy[, "wgs_id"], "-", i)
  metadata_px1 <- rbind(metadata_px1, copy)
}
metadata_px1[, "Sample"] <- metadata_px1[, "wgs_id"]

# Subset the all-Ethiopia IBD to pairs where both samples are Gambella or Mizan
IBD.px1 <- IBD[IBD$sample1 %in% metadata_px1$Sample &
               IBD$sample2 %in% metadata_px1$Sample, ]

edgelist.px1 <- .create.edgelist(IBD.px1)
metadata_px1 <- .reorder.metadata(IBD.px1, metadata_px1)
vertices.px1 <- .create.vertices(IBD.px1, metadata = metadata_px1)

# Per-vertex aesthetics, recomputed after the reorder
px1.labels  <- factor(metadata_px1$px1_haplotype, levels = px1.levels)
px1.colours <- px1.palette[px1.labels]
px1.is.pin  <- metadata_px1$px1_haplotype == "PIN (triple mutant)"
px1.pch     <- pch.map[metadata_px1$px1_shape]
px1.pch[is.na(px1.pch)] <- 21
px1.cex     <- ifelse(px1.is.pin, 2.8, 2)
px1.lwd     <- ifelse(px1.is.pin, 1.5, 0.5)

px1.counts <- table(factor(metadata_px1$px1_haplotype, levels = px1.levels))
px1.keep   <- px1.levels[px1.counts > 0]

district.counts.px1 <- table(factor(metadata_px1$District, levels = px1.districts))

coords.cache.px1 <- vector(mode = "list", length = length(IBD.cutoffs))
names(coords.cache.px1) <- as.character(IBD.cutoffs)

sqrt.n.plots      <- ceiling(sqrt(length(IBD.cutoffs)))
bottom.left.index <- (sqrt.n.plots - 1) * sqrt.n.plots + 1
size              <- 6 * sqrt.n.plots

pdf(file = "IBD_Gambella_Mizan_PX1_PIN.pdf", width = size, height = size)
par(mar = c(1, 11, 1, 1) + 0.1,
    mfrow = c(sqrt.n.plots, sqrt.n.plots),
    cex.main = 1.8)

for (i in seq_along(IBD.cutoffs)) {
  IBD.cutoff <- IBD.cutoffs[i]
  d <- edgelist.px1[edgelist.px1[, "fract_sites_IBD"] >= IBD.cutoff, ]
  IBD.graph <- graph_from_data_frame(d, directed = FALSE, vertices = vertices.px1)

  coords <- coords.cache.px1[[as.character(IBD.cutoff)]]
  if (is.null(coords) || nrow(coords) != vcount(IBD.graph)) {
    coords <- layout_(IBD.graph, nicely())
    coords.cache.px1[[as.character(IBD.cutoff)]] <- coords
  }

  percent.cutoff <- round(IBD.cutoff * 100, digits = 2)

  plot(
    IBD.graph,
    layout       = coords,
    vertex.size  = 6,
    vertex.shape = "none",
    vertex.label = NA,
    main         = paste0("IBD >=", percent.cutoff, "%"),
    rescale      = TRUE
  )

  layout.norm <- norm_coords(coords, xmin = -1, xmax = 1, ymin = -1, ymax = 1)

  # draw PIN carriers last so they sit on top
  draw.order <- order(px1.is.pin)
  points(layout.norm[draw.order, 1], layout.norm[draw.order, 2],
         pch = px1.pch[draw.order], bg = px1.colours[draw.order],
         col = "black", cex = px1.cex[draw.order], lwd = px1.lwd[draw.order])

  if (i == bottom.left.index) {
    legend(
      "bottomleft",
      legend = paste0(px1.keep, " (n=", px1.counts[px1.keep], ")"),
      pch    = 21,
      pt.bg  = px1.palette[px1.keep],
      pt.cex = ifelse(px1.keep == "PIN (triple mutant)", 2, 1.5),
      cex    = 1.2,
      title  = "PfPX1 haplotype",
      inset  = c(-0.25, 0),
      xpd    = NA
    )
    legend(
      "bottomright",
      legend = paste0(px1.districts, " (n=", district.counts.px1[px1.districts], ")"),
      pch    = c(21, 22),
      pt.bg  = "grey50",
      pt.cex = 1.5,
      cex    = 1.2,
      title  = "District",
      inset  = c(-0.1, 0),
      xpd    = NA
    )
  }
}

dev.off()



######## Gambella + Mizan + Uganda: PfPX1 PIN haplotype ########

# Dedicated hmmIBD run over the Gambella, Mizan and Uganda samples, so allele
# frequencies are estimated within this three-population panel rather than
# pooled across all of Ethiopia (unlike the Gambella + Mizan section above).
#
# PfPX1 calls come from px1_genotype_result.tsv, the same source as the section
# above.

IBD.px1u <- read.delim("hmmIBD_Gambella_Mizan_Uganda_maf0.01_out.hmm_fract.txt")

px1u.samples <- unique(c(IBD.px1u$sample1, IBD.px1u$sample2))

# The hmmIBD run defines the sample set, so take metadata for exactly those
# samples rather than re-deriving them with year/fws filters (the Uganda samples
# span 2010-2022, and two of them sit exactly at fws == 0.95).
metadata_px1u <-
  read.csv("Popgen_db_samples_012026_filt_FINAL_fws.csv", fileEncoding = "UTF-8-BOM")
metadata_px1u <- metadata_px1u[metadata_px1u$wgs_id %in% px1u.samples, ]
metadata_px1u$Sample <- metadata_px1u$wgs_id

# Population = the two Ethiopian districts, plus Uganda as a whole (65 samples
# spread thinly over ~8 Ugandan districts, so districts are not resolvable here)
metadata_px1u$px1u_pop <- ifelse(
  metadata_px1u$Country == "Uganda", "Uganda", metadata_px1u$District)

px1u.pops <- c("Gambella", "Mizan", "Uganda")

# --- PfPX1 genotypes: one row per sample, one column per codon ---
# Calls are "WT", "L1222P"/"M1701I"/"D1705N", or NA. read.delim turns the
# numeric headers into X1222 / X1701 / X1705.
px1u <- read.delim("px1_genotype_result.tsv")

px1u.idx <- match(metadata_px1u$wgs_id, px1u$Sample)
metadata_px1u$px1_1222 <- px1u$X1222[px1u.idx]
metadata_px1u$px1_1701 <- px1u$X1701[px1u.idx]
metadata_px1u$px1_1705 <- px1u$X1705[px1u.idx]

px1u.codons <- metadata_px1u[, c("px1_1222", "px1_1701", "px1_1705")]
n.called.u  <- rowSums(!is.na(px1u.codons))
n.mutant.u  <- rowSums(!is.na(px1u.codons) & px1u.codons != "WT")

metadata_px1u$px1_haplotype <- ifelse(
  n.called.u < 3, "Incomplete call",
  ifelse(n.mutant.u == 3, "PIN (triple mutant)",
  ifelse(n.mutant.u == 0, "LMD (wild type)", "Partial (1-2 mutations)")))

# Amino-acid haplotype string (e.g. "LMD", "PIN", "LID"), for the console summary
px1u.wild <- c("L", "M", "D")             # LMD haplotype; mutant is PIN
metadata_px1u$px1_hap_string <- apply(
  mapply(function(call, ref)
    ifelse(is.na(call), "?", ifelse(call == "WT", ref, substr(call, nchar(call), nchar(call)))),
    px1u.codons, px1u.wild),
  1, paste, collapse = "")

px1u.levels <- c("LMD (wild type)", "Partial (1-2 mutations)",
                 "PIN (triple mutant)", "Incomplete call")
px1u.palette <- c(
  "LMD (wild type)"         = "#A6CEE3",
  "Partial (1-2 mutations)" = "#E69F00",
  "PIN (triple mutant)"     = "#D7191C",
  "Incomplete call"         = "#FFFFFF"
)

# Population is encoded by shape
px1u.shape.map <- c("Gambella" = "circle", "Mizan" = "square", "Uganda" = "triangle")
metadata_px1u$px1u_shape <- px1u.shape.map[metadata_px1u$px1u_pop]
metadata_px1u$px1u_shape[is.na(metadata_px1u$px1u_shape)] <- "circle"

px1u.pop.palette <- c("Gambella" = "#E69F00", "Mizan" = "#56B4E9", "Uganda" = "#009E73")

# Build graph objects
edgelist.px1u <- .create.edgelist(IBD.px1u)
metadata_px1u <- .reorder.metadata(IBD.px1u, metadata_px1u)
vertices.px1u <- .create.vertices(IBD.px1u, metadata = metadata_px1u)

# Per-vertex aesthetics, recomputed after the reorder
px1u.labels  <- factor(metadata_px1u$px1_haplotype, levels = px1u.levels)
px1u.colours <- px1u.palette[px1u.labels]
px1u.is.pin  <- metadata_px1u$px1_haplotype == "PIN (triple mutant)"
px1u.pch     <- pch.map[metadata_px1u$px1u_shape]
px1u.pch[is.na(px1u.pch)] <- 21
px1u.cex     <- ifelse(px1u.is.pin, 2.8, 2)
px1u.lwd     <- ifelse(px1u.is.pin, 1.5, 0.5)

px1u.pop.colours <- px1u.pop.palette[metadata_px1u$px1u_pop]
px1u.pop.colours[is.na(px1u.pop.colours)] <- "#FFFFFF"
# second plot: shape carries PIN status instead of population
px1u.pin.pch <- ifelse(px1u.is.pin, 23, 21)

px1u.counts     <- table(factor(metadata_px1u$px1_haplotype, levels = px1u.levels))
px1u.keep       <- px1u.levels[px1u.counts > 0]
px1u.pop.counts <- table(factor(metadata_px1u$px1u_pop, levels = px1u.pops))
px1u.pin.counts <- table(factor(metadata_px1u$px1u_pop, levels = px1u.pops)[px1u.is.pin])

coords.cache.px1u <- vector(mode = "list", length = length(IBD.cutoffs))
names(coords.cache.px1u) <- as.character(IBD.cutoffs)

sqrt.n.plots      <- ceiling(sqrt(length(IBD.cutoffs)))
bottom.left.index <- (sqrt.n.plots - 1) * sqrt.n.plots + 1
size              <- 6 * sqrt.n.plots

.px1u.coords <- function(IBD.graph, IBD.cutoff) {
  coords <- coords.cache.px1u[[as.character(IBD.cutoff)]]
  # only reuse a cached layout if it matches this graph's vertex set
  if (is.null(coords) || nrow(coords) != vcount(IBD.graph)) {
    coords <- layout_(IBD.graph, nicely())
    coords.cache.px1u[[as.character(IBD.cutoff)]] <<- coords
  }
  coords
}


# --- Plot 1: colour = PfPX1 haplotype, shape = population ---
pdf(file = "IBD_Gambella_Mizan_Uganda_PX1_PIN.pdf", width = size, height = size)
par(mar = c(1, 11, 1, 1) + 0.1,
    mfrow = c(sqrt.n.plots, sqrt.n.plots),
    cex.main = 1.8)

for (i in seq_along(IBD.cutoffs)) {
  IBD.cutoff <- IBD.cutoffs[i]
  d <- edgelist.px1u[edgelist.px1u[, "fract_sites_IBD"] >= IBD.cutoff, ]
  IBD.graph <- graph_from_data_frame(d, directed = FALSE, vertices = vertices.px1u)

  coords <- .px1u.coords(IBD.graph, IBD.cutoff)
  percent.cutoff <- round(IBD.cutoff * 100, digits = 2)

  plot(
    IBD.graph,
    layout       = coords,
    vertex.size  = 6,
    vertex.shape = "none",
    vertex.label = NA,
    main         = paste0("IBD >=", percent.cutoff, "%"),
    rescale      = TRUE
  )

  layout.norm <- norm_coords(coords, xmin = -1, xmax = 1, ymin = -1, ymax = 1)

  # draw PIN carriers last so they sit on top
  draw.order <- order(px1u.is.pin)
  points(layout.norm[draw.order, 1], layout.norm[draw.order, 2],
         pch = px1u.pch[draw.order], bg = px1u.colours[draw.order],
         col = "black", cex = px1u.cex[draw.order], lwd = px1u.lwd[draw.order])

  if (i == bottom.left.index) {
    legend(
      "bottomleft",
      legend = paste0(px1u.keep, " (n=", px1u.counts[px1u.keep], ")"),
      pch    = 21,
      pt.bg  = px1u.palette[px1u.keep],
      pt.cex = ifelse(px1u.keep == "PIN (triple mutant)", 2, 1.5),
      cex    = 1.2,
      title  = "PfPX1 haplotype",
      inset  = c(-0.25, 0),
      xpd    = NA
    )
    legend(
      "bottomright",
      legend = paste0(px1u.pops, " (n=", px1u.pop.counts[px1u.pops], ")"),
      pch    = c(21, 22, 24),
      pt.bg  = "grey50",
      pt.cex = 1.5,
      cex    = 1.2,
      title  = "Population",
      inset  = c(-0.1, 0),
      xpd    = NA
    )
  }
}

dev.off()


# --- Plot 2: colour = population, shape = PIN carriage ---
pdf(file = "IBD_Gambella_Mizan_Uganda_PX1_PIN_bypop.pdf", width = size, height = size)
par(mar = c(1, 11, 1, 1) + 0.1,
    mfrow = c(sqrt.n.plots, sqrt.n.plots),
    cex.main = 1.8)

for (i in seq_along(IBD.cutoffs)) {
  IBD.cutoff <- IBD.cutoffs[i]
  d <- edgelist.px1u[edgelist.px1u[, "fract_sites_IBD"] >= IBD.cutoff, ]
  IBD.graph <- graph_from_data_frame(d, directed = FALSE, vertices = vertices.px1u)

  coords <- .px1u.coords(IBD.graph, IBD.cutoff)
  percent.cutoff <- round(IBD.cutoff * 100, digits = 2)

  plot(
    IBD.graph,
    layout       = coords,
    vertex.size  = 6,
    vertex.shape = "none",
    vertex.label = NA,
    main         = paste0("IBD >=", percent.cutoff, "%"),
    rescale      = TRUE
  )

  layout.norm <- norm_coords(coords, xmin = -1, xmax = 1, ymin = -1, ymax = 1)

  draw.order <- order(px1u.is.pin)
  points(layout.norm[draw.order, 1], layout.norm[draw.order, 2],
         pch = px1u.pin.pch[draw.order], bg = px1u.pop.colours[draw.order],
         col = "black", cex = px1u.cex[draw.order], lwd = px1u.lwd[draw.order])

  if (i == bottom.left.index) {
    legend(
      "bottomleft",
      legend = paste0(px1u.pops, " (n=", px1u.pop.counts[px1u.pops], ")"),
      pch    = 21,
      pt.bg  = px1u.pop.palette[px1u.pops],
      pt.cex = 1.5,
      cex    = 1.2,
      title  = "Population",
      inset  = c(-0.25, 0),
      xpd    = NA
    )
    legend(
      "bottomright",
      legend = c(paste0("PIN triple mutant (n=", sum(px1u.is.pin), ")"),
                 paste0("Other / not PIN (n=", sum(!px1u.is.pin), ")")),
      pch    = c(23, 21),
      pt.bg  = "grey50",
      pt.cex = c(2, 1.5),
      cex    = 1.2,
      title  = "PfPX1",
      inset  = c(-0.1, 0),
      xpd    = NA
    )
  }
}

dev.off()


# --- Boxplot: pairwise IBD among PIN carriers vs non-carriers, per population ---
# .build.ibd.by.status / .plot.ibd.boxplot facet on a column called "District",
# so the population label is copied into it here.
metadata_px1u_bp <- metadata_px1u[, c("Sample", "px1u_pop", "px1_haplotype")]
metadata_px1u_bp$District <- metadata_px1u_bp$px1u_pop
metadata_px1u_bp$px1_pin  <- ifelse(
  metadata_px1u_bp$px1_haplotype == "PIN (triple mutant)", "Yes",
  ifelse(metadata_px1u_bp$px1_haplotype == "Incomplete call", NA, "No"))

ibd.df.px1u <- .build.ibd.by.status(IBD.px1u, metadata_px1u_bp, "px1_pin")
.plot.ibd.boxplot(
  ibd.df.px1u,
  "Pairwise IBD: PfPX1 PIN triple mutant vs non-carrier by population",
  "IBD_boxplot_PX1_PIN_Gambella_Mizan_Uganda.pdf",
  metadata_px1u_bp,
  "px1_pin"
)


# --- PfPX1 haplotype summary (printed to console) ---
cat("\nGambella + Mizan + Uganda: PfPX1 haplotype by population\n")
print(table(haplotype  = metadata_px1u$px1_hap_string,
            population = factor(metadata_px1u$px1u_pop, levels = px1u.pops),
            useNA = "ifany"))
print(table(class      = factor(metadata_px1u$px1_haplotype, levels = px1u.levels),
            population = factor(metadata_px1u$px1u_pop, levels = px1u.pops)))

