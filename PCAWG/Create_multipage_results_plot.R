#install.packages(c("pdftools", "magick", "ggplot2", "gridExtra", "qpdf"))
library(pdftools)
library(magick)
library(ggplot2)

# Convert PDF pages to images
image_list = list()
data_dir = '~/dati_Orfeo/scocomello/material_tickTack/PCAWG/results_whole/'
samples=list.files(data_dir)
png_files = lapply(samples[40:60], function(i){
  pdf_file_s <- paste0(data_dir, i, '/plots/plot_timing.png')
  pdf_file_h <- paste0(data_dir, i, '/plots/plot_timing_h.png')
  c(pdf_file_s, pdf_file_h)
}) %>% unlist()

# Function to add a title to a PNG
add_title_to_png <- function(image_path) {
  if (file.exists(image_path)){
  # Read the image
  img <- image_read(image_path)
  title_text = strsplit(image_path, '/')[[1]][7]
  # Convert magick image to raster
  raster_grob <- rasterGrob(as.raster(img), interpolate = TRUE)
  #raster_img + geom_tile(title_text)
  
  # Create ggplot with title overlay
  ggplot() +
    annotation_custom(raster_grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
    ggtitle(title_text) +
    theme_void() +  # Removes unnecessary axes
    theme(plot.title = element_text(size = 6, hjust = 0.5, face = "bold"),
          plot.margin = margin(t = 2, r = 10, b = 5, l = 10))
  }
}

# Apply function to each PNG file
plots <- lapply(seq_along(png_files), function(i) {
  add_title_to_png(png_files[i])
})

# Define output PDF file
output_pdf <- "~/Desktop/multiple_plots.pdf"

# Open PDF device
pdf(output_pdf)  # Adjust size as needed

# Loop through the plot list and print each plot to a new page
for (p in plots) {
  print(p)
}

# Close the PDF device
dev.off()

cat("PDF saved as:", output_pdf)
