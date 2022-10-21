#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(argparser))
suppressPackageStartupMessages(library(highcharter))
suppressPackageStartupMessages(library(htmlwidgets))
suppressPackageStartupMessages(library(ggthemes))
suppressPackageStartupMessages(library(Cairo))


## argparser command line inputs ----

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
  default = NA,
  help = "path to plate reader metadata .xlsx file,
                  default assumes no metadata and just graphs based on Wells"
)

p <- add_argument(p, "--lux",
  short = "-l",
  type = "logical",
  # default = FALSE,
  flag = TRUE,
  help = "Use if running lux (turns on RLU visualization)"
)

p <- add_argument(p, "--bc34",
                  short = "-b",
                  type = "logical",
                  # default = FALSE,
                  flag = TRUE,
                  help = "Use if running Bc34 biosensor (turns on RFI visualization)"
)

p <- add_argument(p, "--output",
  short = "-o",
  default = "output",
  help = "path and prefix for script output files. i.e. 220903_biolog. (where the files go and if they should have anything added to name at beginning)"
)


p <- add_argument(p, "--noGrowth",
                  # default = FALSE,
                  type = "logical",
                  flag = TRUE,
                  help = "Don't make OD graphs (On by default): Include this to turn OFF OD visualization"
)

p <- add_argument(p, "--OD450", 
                  short = "-4",
                  # default = FALSE,
                  type = "logical",
                  flag = TRUE,
                  help = "False by default, assumes OD is \"Absorbance @ 600 OD-600(1)\", turn this on if OD is \"Absorbance @ 450(1)\"."
)


p <- add_argument(p, "--n1to10",
                  # default = FALSE,
                  type = "logical",
                  flag = TRUE,
                  help = "If you using the 10-strains one column per strain layout GK uses, this will add samples via the 1-10 numbering setup (see protocol if unclear)"
)


# Parse the command line arguments
argv <- parse_args(p)

arg_message <- c(
  "\n------------\n",
  "RUNNING SCRIPT \n",
  paste0("INPUT DATA PATH: ", argv$data),
  paste0("METADATA PATH: ", argv$metadata),
  paste0("OUTPUT PATH: ", argv$output),
  paste0("lux=", argv$lux, "  bc34=", argv$bc34, "  OD=", !argv$noGrowth), 
  "\n------------\n"
  )

# print(paste0("input data path is: ", argv$data))
# print(paste0("metadata path is: ", argv$metadata))
# print(paste0("output path is: ", argv$output))
# print(paste0("lux:", argv$lux, " --- bc34:", argv$bc34))

writeLines(arg_message, sep = "\n")


writeLines("input data:")
in_data <- read_csv(argv$data, na = c("", "NA", "N/A"), ) %>%
  filter(!is.na(Label))

metadata_path <- argv$metadata
out_path <- argv$output

## GGPLOT THEME SET AND DEFAULTS
theTheme <- theme_set(ggthemes::theme_clean(base_size = 8))

## functions ----
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
      time_s = hms::round_hms(hms::as_hms(repeatTime * PlateRepeat - repeatTime), 
                              secs = rounding_seconds)) %>%
    mutate(hour = lubridate::hour(time_s)) %>%
    pivot_wider(names_from = measurement, values_from = Result) %>%
    ungroup()
  
  

  if (!is.na(metadata_path)) {
    
    meta_data <- readxl::read_xlsx(metadata_path)
    
    ## If the user used the 10-column layout, attach numbers to the 96-well plate data
    if (argv$n1to10 == TRUE) {
      
      d <- d %>% 
        plater::add_plate(file = "meta_strainNumbers_1to10.csv", 
                          well_ids_column = "Well") %>% 
        mutate(number = factor(number))
      
      meta_data <- meta_data %>% 
        mutate(number = factor(number),
               sample_name = fct_inorder(sample_name))
      }
    
    d <- d %>%
      left_join(meta_data)
  }

  return(d)
}


pt.x_hrScale <- scale_x_time(
  name = "Time (hr)",
  labels = scales::label_time(format = "%H"),
  breaks = scales::breaks_pretty()
)



maxMin <- function(data, readout_col, readout_lbl = NULL) {
    
  d_max <- data %>%     
    filter(!is.na({{ readout_col }})) %>%
    group_by(Well) %>%
    summarize(max = max({{ readout_col }}), min = min({{ readout_col }})) %>%
    pivot_longer(c(max, min), names_to = "max_min", 
                 values_to = rlang::englue("{{ readout_col }}")) %>% 
    mutate(Well = fct_inorder(Well))
  
  p <- d_max %>% 
    ggplot(aes(x = Well, y = {{ readout_col }}, fill = max_min)) +
    geom_col(position = "stack") +
    scale_x_discrete(name = NULL) +
    scale_fill_discrete(name = NULL) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top",
      legend.direction = "horizontal"
    ) + 
    labs(y = readout_lbl)
  
  ggsave(
    filename = paste0(out_path, "_", rlang::englue("{{ readout_col }}"), "minMax_all.pdf"), 
                      p, width = 14, height = 5, device = cairo_pdf)
  
  dname <- paste0(rlang::englue("{{ readout_col }}"), "_maxMin")
  d_max <- list(d_max) %>% setNames(dname)
  
  return(d_max)
}

curve_allWells <- function(data, readout_col, readout_lbl = NULL) {
  
  p <- data %>% 
    ggplot(
      aes(time_s, {{ readout_col }},
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
    labs(y = readout_lbl) +
    theme(legend.position = "right", legend.direction = "vertical") 
  
  
  ggsave(
    filename = paste0(out_path, "_", rlang::englue("{{ readout_col }}"), "curve_allWells.pdf"), 
    p, width = 14, height = 5,
    device = cairo_pdf)
  
  name <- paste0(rlang::englue("{{ readout_col }}"), "_curve_allWells")
  p <- list(p) %>% setNames(name)
  
  return(p)
}

curve_sample <- function(data, readout_col, readout_lbl = NULL) {
  
  p <- data %>% 
    drop_na(sample_name) %>% 
    ggplot(
      aes(time_s, {{ readout_col }},
          color = sample_name,
          fill = sample_name
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
    # scale_color_manual(name = NULL, aesthetics = c("color", "fill")) +
    theme(legend.position = "right", legend.direction = "vertical") + 
    labs(y = readout_lbl)
  
  
  ggsave(
    filename = paste0(out_path, "_", rlang::englue("{{ readout_col }}"), "curve_sample.pdf"), 
    p, width = 14, height = 5, 
    device = cairo_pdf)
  
  name <- paste0(rlang::englue("{{ readout_col }}"), "_curve_sample")
  p <- list(p) %>% setNames(name)
  
  return(p)
}

# curve_hcr_wells <- function(d, readout_col) {
#   
#   # x <- c(
#   #   "Well",
#   #   paste0(readout_col),
#   #   "OD",
#   #   "hour"
#   # )
#   # 
#   # y <- sprintf(
#   #   "{point.%s:.2f}", 
#   #   c("Well",
#   #     readout_col,
#   #     "OD",
#   #     "hour"
#   #   ))
#   
#   # tltip <- tooltip_table(x, y)
#   
#   p <- d %>% hchart("line",
#                     hcaes(
#                       x = hour, 
#                       y = {{ readout_col }},
#                       group = Well
#                     ),
#                     showInLegend = F)
#   # ) %>%
#   #   hc_tooltip(
#   #     useHTML = TRUE,
#   #     headerFormat = "",
#   #     pointFormat = tltip
#   #   )
#   
#   p %>%
#     htmlwidgets::saveWidget(
#       filename = paste0(out_path, "_", 
#                         rlang::englue("{{ readout_col }}"), "curve_interactive.html"), 
#       selfcontained = F
#     )
#   
# }


data_list <- list()
plot_list <- list()


## parse the data and run the script ----
d <- envision_parser(in_data,
  metadata_path = metadata_path,
  rt_min = argv$repeatTime
)
  

## Create OD column. 
## If ABS is 600 (the default), rename column. If ABS is 450, rename column accordingly
d <- d %>% 
  rename(OD = if_else(argv$OD450 == F, 
    "Absorbance @ 600 OD-600", 
    "Absorbance @ 450"))




## If the user sets the data to contain Bc34 or Lux, parse and visualize those based on the expected columns
if (argv$bc34 == TRUE) {
  
  readout_lbl <- "RFI (TurboRFP/Amcyan)"
  
  d <- d %>%
    filter(!is.na(amCyan_Fitnat) | !is.na(turboRFP_Fitnat)) %>%
    mutate(RFI = turboRFP_Fitnat / amCyan_Fitnat) 
  
  d_max <- d %>% 
    maxMin(
      readout_col = RFI, readout_lbl = readout_lbl)
  data_list <- data_list %>% append(d_max)
  
  p <- d %>% curve_allWells(
    readout_col = RFI, readout_lbl = readout_lbl)
  plot_list <- p %>% append(plot_list)
  
  if (!is.na(metadata_path)) {
    
    p <- d %>% 
      curve_sample(
        readout_col = RFI, readout_lbl = readout_lbl)
    plot_list <- p %>% append(plot_list)  
    
  }
  
}



if (argv$lux == TRUE) {
  
  readout_lbl <- "RLU (Lum/OD)"
  
  d <- d %>% 
    rename(Lum = "US LUM 96") %>% 
    mutate(RLU = Lum*5*60/OD)
  
  d_max <- d %>% 
    maxMin(
      readout_col = RLU, readout_lbl = readout_lbl)
  data_list <- data_list %>% append(d_max)
  
  p <- d %>% 
    curve_allWells(
      readout_col = RLU, readout_lbl = readout_lbl)
  plot_list <- p %>% append(plot_list)
  
  if (!is.na(metadata_path)) {
    
    p <- d %>% 
      curve_sample(
        readout_col = RLU, readout_lbl = readout_lbl)
    plot_list <- p %>% append(plot_list)  
    
  }
  
  # d %>% curve_hcr_wells(readout_col = RLU)
  
  
}




## Always plot growth regardless, unless explicitly turned off  
if (argv$noGrowth == FALSE){
  
  readout_lbl <- if_else(argv$OD450 == F, 
                        "OD600",
                        "OD450")

  ## MAX OD GRAPH
  d_max <- d %>%
    filter(!is.na(OD)) %>%
    maxMin(
      readout_col = OD, readout_lbl = readout_lbl)
  
  data_list <- data_list %>% append(d_max)
  
  p <- d %>% 
    curve_allWells(
      readout_col = OD, readout_lbl = readout_lbl)
  
  plot_list <- p %>% append(plot_list)
  
  if (!is.na(metadata_path)) {
    
    p <- d %>% 
      curve_sample(
        readout_col = OD, readout_lbl = readout_lbl)
    plot_list <- p %>% append(plot_list)  
  }
}



if (!is.na(metadata_path)) {

  d_sum <- d %>% group_by(PlateNumber, hour, sample_name)

}

if (is.na(metadata_path)) {

  d_sum <- d %>% group_by(PlateNumber, hour, Well)

}


sum_cols <- c("OD", "RLU", "RFI")
d_sum <- d_sum %>% 
  summarise(
    across(starts_with(sum_cols), 
           list(mean = mean, sd = sd))) %>% 
  mutate(across(starts_with(sum_cols), ~ signif(.x, 3)))
  

data_list <- data_list %>% append(list(summary_data = d_sum))

data_list <- data_list %>% append(setNames(list(d), "data"), after = 0)
data_list %>% writexl::write_xlsx(paste0(out_path, "_data_out.xlsx"))
data_list %>% saveRDS(paste0(out_path, "_data_out.RDS"))

plot_list %>% saveRDS(paste0(out_path, "_plots.RDS"))
