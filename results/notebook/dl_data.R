library(here)

## Download zenodo tar archive and move all files to data directory
## archive described here https://doi.org/10.5281/zenodo.1405579

url_path <- "https://zenodo.org/record/1406586/files/data_2018-08-30.tar.gz?download=1"
output_name <- here("results", "notebook", "data_2018-08-30.tar.gz")

download.file(url_path,
              output_name)

untar(output_name, exdir = here("results", "notebook"))
dir.create(here("tmp_data"))

file.rename(here("results", "notebook", "data", 
                 "results_rds",
                 "2018-05-16_mouse_human",
                 "processed_data.rds"),
            here("results", 
                 "2018-05-16_mouse_human",
                 "processed_data.rds"))

file.rename(here("results", "notebook", "data", 
                 "results_rds",
                 "2018-05-15_pbmc",
                 "processed_data.rds"),
            here("results", 
                 "2018-05-15_pbmc",
                 "processed_data.rds"))

unlink(here("results",
            "notebook",
            "results_rds"),
       force = T,
       recursive = T)       

file.rename(here("results", "notebook", "data"), 
            here("data"))

unlink(output_name)
