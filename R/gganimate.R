library(tidyverse)
library(gganimate)
library(ggrepel)
library(magick)
library(cowplot)
library(here)

project_dir <- here::here()
data_dir <- file.path(project_dir, "data")
results_dir <- file.path(project_dir, "results")
color_palette <- c("#0072B2", 
                   "#D55E00")

cells <- list(
  mouse_human_cell_pulldown = c("GACGTTAGTGCCTGTG",
                                "CTGATCCCATGACGGA"))

libs <- c(
  "original_10x",
  "mouse_human_cell_pulldown")

lib_data_dir <- file.path(data_dir, 
                          "lna_cell",
                          "mh_mix")

## original library to compare against
reflib <- "original_10x"
resampled_libs <- c("mouse_human_cell_pulldown")

## reference resampled lib for resampled vs control plots
resampled_lib <- "mouse_human_cell_pulldown"

## pretty name for libraries
lib_names = c(
  original_10x = "Original Library",
  mouse_human_cell_pulldown = "Resampled Library"
)

sc_objs <- readRDS(file.path(results_dir, 
                             "2018-05-16_mouse_human", 
                             "processed_data.rds"))

sc_metadat <- map(sc_objs, ~.x$meta_dat) %>% 
  bind_rows(.id = "library") %>% 
  mutate(library = lib_names[library],
         library = factor(library, levels = lib_names)) %>% 
  arrange(resampled)

plt1 <- ggplot(sc_metadat, aes(total_umis, cDNA_duplication)) +
  geom_point(aes(color = resampled), size = 1) +
  geom_text_repel(data = filter(sc_metadat, resampled),
                  aes(label = ".",
                      color = resampled),
                  size = 0,
                  force = 10,
                  min.segment.length = 0,
                  point.padding = 1.0,
                  seed = 42,
                  arrow = arrow(length = unit(0.3,
                                              "line"), 
                                angle = 35,
                                type = "open", 
                                ends = "last")
  ) +
  labs(x = "# of UMIs",
       y = "Sequencing saturation",
       title = "{next_state}") + 
  scale_x_log10(labels = scales::comma) +
  scale_color_manual(values = color_palette) +
  theme_cowplot(font_size = 16, line_size = 1)  +
  theme(legend.position = "none") +
  transition_states(
    library,
    transition_length = 2,
    state_length = 1
  ) +
  ease_aes('linear')

plt2 <- ggplot(sc_metadat,
               aes(human_umis, 
                   mouse_umis, 
                   colour = resampled)) +
  geom_point(size = 0.5) + 
  geom_text_repel(data = filter(sc_metadat, resampled),
                  aes(label = ".",
                      color = resampled),
                  size = 0,
                  force = 10,
                  min.segment.length = 0,
                  point.padding = 1.15,
                  arrow = arrow(length = unit(0.3,
                                              "line"), 
                                angle = 35,
                                type = "open", 
                                ends = "last"
                  ),
                  seed = 2
  ) + 
  labs(x = "Human UMIs",
       y = "Mouse UMIs",
       title = "{closest_state}") +
  scale_x_continuous(labels = scales::comma) + 
  scale_y_continuous(labels = scales::comma) +
  scale_colour_manual(name = "resampled:",
                      values = color_palette) +
  theme_cowplot(font_size = 16, line_size = 1) +
  theme(legend.position = "none",
        legend.text = element_text(size = 12)) +
  transition_states(
    library,
    transition_length = 2,
    state_length = 1
  ) +
  ease_aes('linear')

plt2

### fig 2c
cells <- list(
  mkcell_pulldown = c("TGCGCAGCAGGTCGTC",
                      "ACTTGTTAGGACCACA",
                      "CCATTCGTCCCTGACT",
                      "TGTCCCAGTAAACACA"))

libs <- c(
  "kirkpatrick",
  "mkcell_pulldown")


lib_data_dir <- file.path(data_dir, 
                          "lna_cell",
                          "pbmc_expt")

## original library to compare against
reflib <- "kirkpatrick"
## all resampled libs to plot
resampled_libs <- "mkcell_pulldown"
## reference resampled lib for resampled vs control plots
resampled_lib <- "mkcell_pulldown"

## pretty name for libraries
lib_names = c(
  kirkpatrick = "Original Library",
  mkcell_pulldown = "Resampled Library"
)

sc_objs <- readRDS(file.path(results_dir, "2018-05-15_pbmc", "processed_data.rds"))

sc_metadat <- map(sc_objs, ~.x$meta_dat) %>% 
  bind_rows(.id = "library") %>% 
  mutate(library = lib_names[library],
         library = factor(library, levels = lib_names)) %>% 
  arrange(resampled)


plt3 <- ggplot(sc_metadat, aes(total_umis, cDNA_duplication)) +
  geom_point(aes(color = resampled), size = 0.5) +
  scale_x_continuous(labels = scales::comma) +
  scale_color_manual(values = color_palette) +
  geom_text_repel(data = filter(sc_metadat, resampled),
                  aes(label = ".",
                      color = resampled),
                  size = 0,
                  force = 10,
                  min.segment.length = 0,
                  point.padding = 1.0,
                  arrow = arrow(length = unit(0.3,
                                              "line"), 
                                angle = 35,
                                type = "open", 
                                ends = "last")
  ) +
  theme_cowplot(font_size = 16, line_size = 1)  +
  theme(legend.position = "none",
        plot.margin = unit(c(5.5, 20.5, 5.5, 5.5), 
                           "points")) +
  labs(x = "# of UMIs",
       y = "Sequencing saturation",
         title = "{closest_state}") + 
  transition_states(
    library,
    transition_length = 2,
    state_length = 1
  ) +
  ease_aes('linear')

plt3


# make multiple gifs into one
# see https://github.com/thomasp85/gganimate/wiki/Animation-Composition

a_plt <- animate(plt1)
b_plt <- animate(plt2)
c_plt <- animate(plt3)

a_gif <- image_read(a_plt)
b_gif <- image_read(b_plt)
c_gif <- image_read(c_plt)

new_gif <- image_append(c(b_gif[1], 
                          a_gif[1],
                          c_gif[1]))

for(i in 2:100){
  combined <- image_append(c(b_gif[i], 
                             a_gif[i],
                             c_gif[i]))
  new_gif <- c(new_gif, combined)
}

new_gif
