ggplot2::theme_set(theme_bw())
suppressWarnings(
theme_update(
  axis.text.x = element_text(size = 16)
, axis.text.y = element_text(size = 16)
, axis.title.x = element_text(size = 16)
, axis.title.y = element_text(size = 16)
, legend.title = element_text(size = 16)
, legend.text  = element_text(size = 14)
, panel.grid.major = element_blank()
, panel.grid.minor = element_blank()
, strip.background = element_blank()
, panel.margin = unit(0, "lines")
, legend.key.size = unit(.55, "cm")
, legend.key = element_rect(fill = "white")
, panel.margin.y = unit(0.5, "lines")
, panel.border = element_rect(
    colour = "black", fill = NA, size = 1)
, strip.text.x = element_text(
    size = 16, colour = "black", face = "bold")
  )
)


