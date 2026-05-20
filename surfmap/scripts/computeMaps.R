rm(list=ls())

###################  LIBRARIES   ##########################
packages <- c("optparse")

installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
suppressMessages(invisible(lapply(packages, library, character.only = TRUE)))
options(warn=-1)

################## FUNCTIONS ##########################

image.nan.better <- function(z, zlim, col,
                             na.color='gray',
                             outside.below.color='black',
                             outside.above.color='blue', ...) {
  zstep <- (zlim[2] - zlim[1]) / length(col)
  newz.below.outside <- zlim[1] - 2 * zstep
  newz.above.outside <- zlim[2] + zstep
  newz.na <- zlim[2] + 2 * zstep

  z[which(z < zlim[1])] <- newz.below.outside
  z[which(z > zlim[2])] <- newz.above.outside
  z[which(is.na(z))] <- newz.na

  zlim[1] <- zlim[1] - 2 * zstep
  zlim[2] <- zlim[2] + 2 * zstep

  col <- c(outside.below.color, col[1], col, outside.above.color, na.color)
  image(z = z, zlim = zlim, col = col, ...)
}

image.scale <- function(z, zlim, col, scalename, breaks,
                        horiz = TRUE, ylim = NULL, xlim = NULL, ...) {
  if (missing(breaks) && !missing(zlim)) {
    breaks <- seq(zlim[1], zlim[2], length.out = (length(col) + 1))
  }
  if (missing(breaks) && missing(zlim)) {
    zlim <- range(z, na.rm = TRUE)
    zlim[2] <- zlim[2] + (zlim[2] - zlim[1]) * (1E-3)
    zlim[1] <- zlim[1] - (zlim[2] - zlim[1]) * (1E-3)
    breaks <- seq(zlim[1], zlim[2], length.out = (length(col) + 1))
  }

  poly <- vector(mode = "list", length(col))
  for (i in seq_along(poly)) {
    poly[[i]] <- c(breaks[i], breaks[i + 1], breaks[i + 1], breaks[i])
  }

  xaxt <- ifelse(horiz, "s", "n")
  yaxt <- ifelse(horiz, "n", "s")

  if (horiz) { YLIM <- c(0, 1); XLIM <- range(breaks) }
  if (!horiz) { YLIM <- range(breaks); XLIM <- c(0, 1) }

  if (is.null(xlim)) xlim <- XLIM
  if (is.null(ylim)) ylim <- YLIM

  plot(1, 1, t = "n", ylim = ylim, xlim = xlim, xaxt = xaxt, yaxt = yaxt,
       xaxs = "i", yaxs = "i", xlab = scalename, ylab = "", ...)

  for (i in seq_along(poly)) {
    if (horiz) polygon(poly[[i]], c(0, 0, 1, 1), col = col[i], border = NA)
    if (!horiz) polygon(c(0, 0, 1, 1), poly[[i]], col = col[i], border = NA)
  }
}

################## MAIN ##########################

option_list <- list(
  make_option(c("-i", "--input"), type="character", default=".", help="path to matrix files or matrix file", metavar="character"),
  make_option(c("-p", "--pdb"), type="character", default=".", help="pdb id (used in title)", metavar="character"),
  make_option(c("-c", "--coord"), type="character", default=NA, help="File containing (phi, theta) coordinates to map", metavar="character"),
  make_option(c("-l", "--reslist"), type="character", default=NA, help="File containing coordinates of residues to map", metavar="character"),
  make_option(c("-s", "--cellsize"), type="integer", default=5, help="grid cell size (must divide 180)", metavar="integer"),
  make_option(c("-P", "--projection"), type="character", default="sinusoidal", help="projection type", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=".", help="output directory", metavar="character"),
  make_option(c("--suffix"), type="character", default="_smoothed_matrix.txt", help="suffix removed from input to create output basename", metavar="character"),
  make_option(c("--png"), action="store_true", default=FALSE, help="output png instead of pdf"),
  make_option(c("--margin_scale"), type="double", default=1.0, help="scale factor applied to par(mar=...). Smaller => less whitespace", metavar="double"),
  make_option(c("--no_scale_bar"), action="store_true", default=FALSE, help="disable the color scale bar panel on the right"),
  make_option(c("--electrostatics"), action="store_true", default=FALSE, help="use electrostatics scale"),
  make_option(c("--stickiness"), action="store_true", default=FALSE, help="use stickiness scale"),
  make_option(c("--kyte_doolittle"), action="store_true", default=FALSE, help="use Kyte-Doolittle scale"),
  make_option(c("--wimley_white"), action="store_true", default=FALSE, help="use Wimley-White scale"),
  make_option(c("--circular_variance"), action="store_true", default=FALSE, help="use circular variance scale"),
  make_option(c("--bfactor"), action="store_true", default=FALSE, help="use bfactor scale"),
  make_option(c("--discrete"), action="store_true", default=FALSE, help="use discrete scale"),
  make_option(c("--elec_max_value"), type="double", default=NULL, help="max abs value for electrostatics scale", metavar="double"),
  make_option(c("--bfactor_min_value"), type="double", default=NULL, help="min value for bfactor scale", metavar="double"),
  make_option(c("--bfactor_max_value"), type="double", default=NULL, help="max value for bfactor scale", metavar="double")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

margin_scale <- suppressWarnings(as.numeric(opt$margin_scale))
if (length(margin_scale) != 1 || is.na(margin_scale) || margin_scale <= 0) {
  margin_scale <- 1.0
}

show_scale_bar <- !(opt$no_scale_bar == TRUE)
TIGHT_MARGIN_THRESHOLD <- 0.01
tight_png <- isTRUE(opt$png) && !show_scale_bar && (margin_scale <= TIGHT_MARGIN_THRESHOLD)

if (file_test("-f", opt$input)) {
  files <- c(opt$input)
} else if (file_test("-d", opt$input)) {
  files <- list.files(opt$input, pattern = "\\matrix.txt$", full.names = TRUE)
} else {
  cat("Error from computeMaps.R: -i is neither a file nor a directory\n")
  quit(status = 1)
}

width <- as.integer(opt$cellsize)
if (180 %% width != 0) {
  cat("Error: cellsize must divide 180\n")
  quit(status = 1)
}
stepabs <- 360 / width
stepord <- 180 / width
asp_map <- stepord / stepabs

outdir <- file.path(opt$outdir, "maps")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

for (f in files) {
  name_prefix <- gsub(opt$suffix, "", basename(f))
  pdf_filename <- file.path(outdir, paste0(name_prefix, "_map.pdf"))

  if (opt$png) {
    if (tight_png) {
      target_w_px <- 2100
      px_per_cell <- max(1, floor(target_w_px / stepabs))
      w_px <- stepabs * px_per_cell
      h_px <- stepord * px_per_cell
      png(gsub("\\.pdf$", ".png", pdf_filename), width = w_px, height = h_px, units = "px")
    } else {
      png(gsub("\\.pdf$", ".png", pdf_filename), res = 300, width = 17.78, height = 17.78, units = "cm")
    }
  } else {
    pdf(pdf_filename)
  }

  par(oma = c(0, 0, 0, 0))
  if (tight_png) {
    par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
  }

  data_matrix <- read.table(f, fill = TRUE, header = TRUE, sep = "\t")
  data_matrix[is.na(data_matrix)] <- 0
  val_matrix <- matrix(as.double(data_matrix[, 3]), ncol = stepabs, nrow = stepord, byrow = FALSE)

  proj <- which(val_matrix != Inf)
  minval <- min(val_matrix[proj])
  maxval <- max(val_matrix[proj])

  if (opt$electrostatics) {
    if (is.null(opt$elec_max_value)) {
      minval <- min(val_matrix[proj])
      maxval <- max(val_matrix[proj])
    } else {
      minval <- -abs(opt$elec_max_value)
      maxval <-  abs(opt$elec_max_value)
    }
    if (abs(minval) > abs(maxval)) {
      maxval <- abs(minval)
    } else {
      minval <- -abs(maxval)
    }
    rangev <- abs(minval - maxval)
    scale_main <- paste0("electrostatic\npotential")
    main_title <- paste0("electrostatic potential map\n", opt$pdb)
    colors <- c(seq(minval, minval + rangev * 1/3, length=334),
                seq(minval + rangev * 1/3, minval + rangev * 2/3, length=333),
                seq(minval + rangev * 2/3, maxval, length=334))
    scale_at <- c(minval, minval + rangev/6, minval + 2*rangev/6, minval + 3*rangev/6,
                  minval + 4*rangev/6, minval + 5*rangev/6, maxval)
    colorScale <- colorRampPalette(c("red", "white", "blue"))(1000)

  } else if (opt$kyte_doolittle) {
    minval <- -4.5; maxval <- 4.5
    rangev <- abs(minval - maxval)
    scale_main <- "hydrophobicity\nKyte-Doolittle"
    main_title <- paste0("Kyte-Doolittle hydrophobicity map\n", opt$pdb)
    colors <- c(seq(minval, minval + rangev * 1/3, length=334),
                seq(minval + rangev * 1/3, minval + rangev * 2/3, length=333),
                seq(minval + rangev * 2/3, maxval, length=334))
    scale_at <- c(-4.5,-3,-1.5,0,1.5,3,4.5)
    colorScale <- colorRampPalette(c("cadetblue", "cadetblue3", "#faf4e0", "orange3", "sienna4"))(1000)

  } else if (opt$stickiness) {
    minval <- -1.273; maxval <- 1.273
    rangev <- abs(minval - maxval)
    scale_main <- "stickiness"
    main_title <- paste0("stickiness map\n", opt$pdb)
    colors <- c(seq(minval, minval + rangev * 1/3, length=334),
                seq(minval + rangev * 1/3, minval + rangev * 2/3, length=333),
                seq(minval + rangev * 2/3, maxval, length=334))
    scale_at <- c(-1.273,-0.85,-0.43,0,0.43,0.85,1.273)
    colorScale <- colorRampPalette(c("royalblue3", "white", "darkgreen"))(1000)

  } else if (opt$circular_variance) {
    minval <- 0; maxval <- 1
    rangev <- 1
    scale_main <- "circular\nvariance"
    main_title <- paste0("circular variance map\n", opt$pdb)
    colors <- c(seq(0, 1, length=1000))
    scale_at <- c(0,1/6,2/6,3/6,4/6,5/6,1)
    colorScale <- colorRampPalette(c("black", "white", "blue"))(1000)

  } else if (opt$bfactor) {
    minval <- min(val_matrix[proj]); maxval <- max(val_matrix[proj])
    if (!is.null(opt$bfactor_min_value)) minval <- opt$bfactor_min_value
    if (!is.null(opt$bfactor_max_value)) maxval <- opt$bfactor_max_value
    rangev <- abs(minval - maxval)
    scale_main <- "b-factor"
    main_title <- paste0("b-factor map\n", opt$pdb)
    colors <- c(seq(minval, minval + rangev * 1/3, length=334),
                seq(minval + rangev * 1/3, minval + rangev * 2/3, length=333),
                seq(minval + rangev * 2/3, maxval, length=334))
    scale_at <- c(minval, minval + rangev/6, minval + 2*rangev/6, minval + 3*rangev/6,
                  minval + 4*rangev/6, minval + 5*rangev/6, maxval)
    colorScale <- colorRampPalette(c("chocolate3", "white", "darkblue"))(1000)

  } else {
    rangev <- abs(minval - maxval)
    scale_main <- "value"
    main_title <- paste0("map\n", opt$pdb)
    colors <- c(seq(minval, maxval, length=1000))
    scale_at <- c(minval, minval + rangev/6, minval + 2*rangev/6, minval + 3*rangev/6,
                  minval + 4*rangev/6, minval + 5*rangev/6, maxval)
    colorScale <- colorRampPalette(c("blue", "white", "red"))(1000)
  }

  if (show_scale_bar) {
    layout(matrix(c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))
  } else {
    layout(matrix(1, nrow=1, ncol=1))
  }

  if (!tight_png) {
    par(mar = c(14.4, 5, 9.3, 1.2) * margin_scale)
  }

  if (tight_png) {
    labx <- ""
    laby <- ""
  } else {
    labx <- expression(paste(phi, " sin(", theta, ")"))
    laby <- expression(paste("90 - ", theta))
  }

  # Build target Matrix for tracing MAb tags perimeters
  tag_matrix <- matrix(NA, nrow=stepord, ncol=stepabs)
  if (!is.na(opt$reslist) && file.exists(opt$reslist)) {
      res_data <- tryCatch(read.table(opt$reslist, fill=TRUE, header=FALSE, stringsAsFactors=FALSE), error=function(e) NULL)
      if (!is.null(res_data) && ncol(res_data) >= 4 && any(grepl("CDR", res_data[,4], ignore.case=TRUE))) {
          
          target_residues <- paste(res_data[,3], res_data[,2], res_data[,1], sep="_")
          target_tags <- res_data[,4]

          tag_vec <- rep(NA, nrow(data_matrix))

          for(i in 1:nrow(data_matrix)) {
              if (!is.na(data_matrix$residues[i]) && data_matrix$residues[i] != "") {
                  cell_res <- trimws(unlist(strsplit(as.character(data_matrix$residues[i]), ",")))
                  matches <- cell_res[cell_res %in% target_residues]
                  if (length(matches) > 0) {
                      idx <- match(matches[1], target_residues)
                      tag_vec[i] <- target_tags[idx]
                  }
              }
          }
          tag_matrix <- matrix(tag_vec, nrow=stepord, ncol=stepabs, byrow=FALSE)
      }
  }
  t_tag_matrix <- t(tag_matrix)

  image.nan.better(t(val_matrix), col=colorScale, zlim=c(minval, maxval), outside.below.color='white', outside.above.color='gray90', na.color='white', frame.plot=!tight_png, axes=FALSE, xlab=labx, ylab=laby, asp=asp_map, cex.lab=1.5)

  # Draw perimeter borders for MAb tags
  if (exists("t_tag_matrix") && any(!is.na(t_tag_matrix))) {
      dx <- 1 / (stepabs - 1)
      dy <- 1 / (stepord - 1)

      get_color <- function(tag) {
          tag <- toupper(tag)
          if(tag == "CDR3H") return("darkred")
          if(tag == "CDR2H") return("red")
          if(tag == "CDR1H") return("lightcoral")
          if(tag == "CDR3L") return("darkblue")
          if(tag == "CDR2L") return("blue")
          if(tag == "CDR1L") return("deepskyblue")
          return("black")
      }

      get_lwd <- function(tag) {
          tag <- toupper(tag)
          if(tag == "CDR3H" || tag == "CDR3L") return(3.5)
          if(tag == "CDR2H" || tag == "CDR2L") return(2.5)
          if(tag == "CDR1H" || tag == "CDR1L") return(1.5)
          return(2.5)
      }

      for(i in 1:stepabs) {
          for(j in 1:stepord) {
              tag <- t_tag_matrix[i, j]
              if(!is.na(tag)) {
                  col <- get_color(tag)
                  lwd_val <- get_lwd(tag)
                  
                  x_center <- (i - 1) * dx
                  y_center <- (j - 1) * dy

                  x_left <- x_center - dx/2
                  x_right <- x_center + dx/2
                  y_bottom <- y_center - dy/2
                  y_top <- y_center + dy/2

                  if(i == 1 || is.na(t_tag_matrix[i-1, j]) || t_tag_matrix[i-1, j] != tag) {
                      segments(x_left, y_bottom, x_left, y_top, col=col, lwd=lwd_val, lty=3)
                  }
                  if(i == stepabs || is.na(t_tag_matrix[i+1, j]) || t_tag_matrix[i+1, j] != tag) {
                      segments(x_right, y_bottom, x_right, y_top, col=col, lwd=lwd_val, lty=3)
                  }
                  if(j == 1 || is.na(t_tag_matrix[i, j-1]) || t_tag_matrix[i, j-1] != tag) {
                      segments(x_left, y_bottom, x_right, y_bottom, col=col, lwd=lwd_val, lty=3)
                  }
                  if(j == stepord || is.na(t_tag_matrix[i, j+1]) || t_tag_matrix[i, j+1] != tag) {
                      segments(x_left, y_top, x_right, y_top, col=col, lwd=lwd_val, lty=3)
                  }
              }
          }
      }
  }

  if (!tight_png) {
    axis(1, at=c(0,0.25,0.5,0.75,1), labels=c(-180,-90,0,90,180), cex.axis=1.2)
    axis(2, at=c(0,0.25,0.5,0.75,1), labels=c(-90,-45,0,45,90), cex.axis=1.2, las=2)
    title(main = main_title, line = 1.5)
  }

  if (show_scale_bar) {
    par(mar = c(15.5, 1.6, 10.5, 4.5) * margin_scale)
    image.scale(t(val_matrix), col=colorScale, breaks=colors, scalename=scale_main, horiz=FALSE, yaxt="n")
    axis(4, at=scale_at, las=2, cex.axis=0.8, labels=round(scale_at, digits=2))
  }

  dev.off()
}

quit(status=0)
