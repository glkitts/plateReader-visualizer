#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(argparser))
suppressPackageStartupMessages(library(highcharter))
suppressPackageStartupMessages(library(htmlwidgets))
suppressPackageStartupMessages(library(ggthemes))
suppressPackageStartupMessages(library(Cairo))


## Command line inputs ----

# Create a parser
p <- argparser::arg_parser("Plate Reader Data Analyzer-Visualizer")

# Add command line arguments
p <- add_argument(p, "data",
  type = "character",
  help = "path to plate reader input csv file (listA format) with no plate info included (Well data only)"
)

p <- add_argument(p, "--repeatTime",
  short = "-r",
  default = 15,
  help = "time between measurement of same well (in minutes)",
  type = "numeric"
)

p <- add_argument(p, "--metadata",
  short = "-m",
  help = "path to plate reader metadata .xlsx file,
                  default assumes no metadata and just graphs based on Wells"
)

p <- add_argument(p, "--n1to10",
  default = FALSE,
  type = "logical",
  help = "If you using the 10-strains one column per strain layout GK uses, this will add samples via the 1-10 numbering setup (see protocol if unclear)"
)

p <- add_argument(p, "--lux",
  short = "-l",
  default = FALSE,
  type = "logical",
  help = "Use if running lux (turns on RLU visualization)"
)

p <- add_argument(p, "--bc34",
                  short = "-b",
                  default = FALSE,
                  type = "logical",
                  help = "Use if running Bc34 biosensor (turns on RFI visualization)"
)

p <- add_argument(p, "--noGrowth",
                  default = TRUE,
                  type = "logical",
                  help = "Makes OD graphs by default, use this to turn off growth visualization"
)



p <- add_argument(p, "--output",
  short = "-o",
  default = "output",
  help = "path and prefix for script output files. i.e. 220903_biolog. (where the files go and if they should have anything added to name at beginning)"
)


# Parse the command line arguments
argv <- parse_args(p)




in_data <- read_csv(argv$data, na = c("", "NA", "N/A")) %>%
  filter(!is.na(Label))

metadata_path <- argv$metadata
out_path <- argv$output

## GGPLOT THEME SET AND DEFAULTS
theTheme <- theme_set(ggthemes::theme_clean(base_size = 8))

## parse the plate reader data ----
envision_parser <- function(in_data,
                                 metadata_path = NA,
                                 parse_meta = T,
                                 rt_min = 15,
                                 rounding_seconds = 1800) {
  
  
  repeatTime <- hms::hms(minutes = rt_min)


  d <- in_data %>%
    filter(Label != "N/A") %>%
    group_by(PlateRepeat, Well) %>%
    select(PlateNumber, PlateRepeat, Well, Label, Result) %>%
    separate(Label, into = c("measurement", "group"), sep = "\\(") %>%
    mutate(
      group = str_remove(group, "\\)"),
      time = repeatTime * PlateRepeat - repeatTime
    ) %>%
    mutate(hour = (PlateRepeat * rt_min - rt_min) / 60) %>%
    pivot_wider(names_from = measurement, values_from = Result) %>%
    ungroup()

  if (!is.na(metadata_path)) {
    meta_data <- readxl::read_xlsx(metadata_path) %>%
      mutate(PlateNumber = factor(PlateNumber))
    d <- d %>%
      mutate(PlateNumber = factor(PlateNumber)) %>%
      left_join(meta_data)
  }

  return(d)
}


pt.x_hrScale <- scale_x_time(
  name = "Time (hr)",
  labels = scales::label_time(format = "%H"),
  breaks = scales::breaks_pretty()
)

## parse the data and run the script ----

d <- envision_parser(in_data,
  metadata_path = metadata_path,
  rt_min = argv$repeatTime
) %>%
  rename(OD = "Absorbance @ 450")


d_RFI <- d %>%
  filter(!is.na(amCyan_Fitnat) | !is.na(turboRFP_Fitnat)) %>%
  mutate(RFI = turboRFP_Fitnat / amCyan_Fitnat)

d_RFI %>% colnames() # readr::write_csv("./dRFI_csv")


# MAX RFI BAR GRAPH
p.maxRFI <- d_RFI %>%
  group_by(Well) %>%
  summarize(max = max(RFI), min = min(RFI)) %>%
  pivot_longer(c(max, min), names_to = "max_min", values_to = "RFI") %>%
  # arrange(desc(Well)) %>%
  mutate(Well = fct_inorder(Well)) %>%
  ggplot(aes(x = Well, y = RFI, fill = max_min)) +
  geom_col(position = "stack") +
  scale_x_discrete(name = NULL) +
  scale_y_continuous(expand = c(0.001, 0.001)) +
  scale_fill_discrete(name = NULL) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    legend.direction = "horizontal"
  )


# MAX RFI BAR GRAPH OUPUT SETTINGS
ggsave(filename = paste0(out_path, "_RFIminMax_all.pdf"), p.maxRFI, width = 14, height = 5)




# MAX OD GRAPH
# d <- d %>% rename(OD = "Absorbance @ 450")

p.maxOD <- d %>%
  filter(!is.na(OD)) %>%
  group_by(Well) %>%
  summarize(max = max(OD), min = min(OD)) %>%
  pivot_longer(c(max, min), names_to = "max_min", values_to = "OD450") %>%
  mutate(Well = fct_inorder(Well)) %>%
  ggplot(aes(x = Well, y = OD450, fill = max_min)) +
  geom_col(position = "stack") +
  scale_x_discrete(name = NULL) +
  scale_y_continuous(expand = c(0.001, 0.001)) +
  scale_fill_discrete(name = NULL) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    legend.direction = "horizontal"
  )

# MAX OD GRAPH OUPUT SETTINGS
ggsave(filename = paste0(out_path, "_ODminMax_all.pdf"), p.maxOD, width = 14, height = 5)



## HighCharter Interactive Plots ----

x <- c(
  "Well",
  "RFI",
  "OD",
  "hour"
)

y <- sprintf(
  "{point.%s:.2f}",
  c(
    "Well",
    "RFI",
    "OD",
    "hour"
  )
)
tltip <- tooltip_table(x, y)

### RFI curve for all wells ----
p.hcr_RFIcurve_Well <- d_RFI %>%
  hchart("line",
    hcaes(
      x = hour, y = RFI,
      group = Well
    ),
    showInLegend = F
  ) %>%
  hc_tooltip(
    useHTML = TRUE,
    headerFormat = "",
    pointFormat = tltip
  )

p.hcr_RFIcurve_Well %>%
  htmlwidgets::saveWidget(
    file = paste0("./", out_path, "_RFI_perWell_interactive.html"),
    selfcontained = F
  )

### OD curve for all wells ----
p.hcr_ODcurve_Well <- d_RFI %>%
  hchart("line",
    hcaes(
      x = hour, y = OD,
      group = Well
    ),
    showInLegend = F
  ) %>%
  hc_tooltip(
    useHTML = TRUE,
    headerFormat = "",
    pointFormat = tltip
  )

p.hcr_ODcurve_Well %>%
  htmlwidgets::saveWidget(
    file = paste0("./", out_path, "_OD_perWell_interactive.html"),
    selfcontained = F
  )




# RFI ALL WELLS CURVE
p.RFIcurve_Wells <- d_RFI %>%
  ggplot(aes(time, RFI,
    color = Well,
    fill = Well
  )) +
  stat_summary(
    fun.data = "mean_sdl",
    fun.args = list(mult = 1),
    geom = "ribbon",
    color = NA,
    alpha = 0.2
  ) +
  pt.x_hrScale +
  stat_summary(
    fun = "mean",
    geom = "line"
  ) +
  theme(legend.position = "right", legend.direction = "vertical")

# RFI ALL WELLS CURVE OUTPUT SETTINGS
ggsave(paste0("./", out_path, "_RFIcurve_allWells.pdf"), plot = p.RFIcurve_Wells, width = 17, height = 11)


if (!is.na(metadata_path)) {
  d_RFI_samples <- d_RFI %>%
    filter(!is.na(sample_name)) %>%
    mutate(
      sample_name = fct_inorder(sample_name),
      color = fct_inorder(color)
    )

  p.RFIcurve_sample <- d_RFI_samples %>%
    ggplot(aes(time, RFI,
      color = color,
      fill = color
    )) +
    stat_summary(
      fun.data = "mean_sdl",
      fun.args = list(mult = 1),
      geom = "ribbon",
      color = NA,
      alpha = 0.2
    ) +
    stat_summary(
      fun = "mean",
      geom = "line"
    ) +
    scale_color_identity(
      name = NULL,
      labels = levels(d_RFI_samples$sample_name),
      breaks = levels(d_RFI_samples$color),
      aesthetics = c("fill", "color"),
      guide = "legend"
    ) +
    pt.x_hrScale +
    labs(y = "RFI (TurboRFP/Amcyan)")

  ggsave(paste0("./", out_path, "_RFIcurve_samples.pdf"), plot = p.RFIcurve_sample, width = 7, height = 4, device = cairo_pdf)

  p.ODcurve_sample <- d_RFI_samples %>%
    ggplot(aes(time, OD,
      color = color,
      fill = color
    )) +
    stat_summary(
      fun.data = "mean_sdl",
      fun.args = list(mult = 1),
      geom = "ribbon",
      color = NA,
      alpha = 0.2
    ) +
    stat_summary(
      fun = "mean",
      geom = "line"
    ) +
    scale_color_identity(
      name = NULL,
      labels = levels(d_RFI_samples$sample_name),
      breaks = levels(d_RFI_samples$color),
      aesthetics = c("fill", "color"),
      guide = "legend"
    ) +
    pt.x_hrScale +
    scale_y_continuous(trans = "log2")


  ggsave(paste0("./", out_path, "_ODcurve_samples.pdf"), plot = p.ODcurve_sample, width = 7, height = 4, device = cairo_pdf)
}
