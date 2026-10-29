# Correct only degenerate lobster pools; preserve every non-singleton cached draw.
# Rscript code/R/correct_lobster_singleton_cache.R BACKUP_DIRECTORY
local({
  suppressPackageStartupMessages({library(dplyr); library(here)})
  source(here::here("code/R/01_utils.R"),local=TRUE)
  args <- commandArgs(trailingOnly=TRUE);stopifnot(length(args)==1L)
  dir.create(args[1],recursive=TRUE,showWarnings=FALSE)
  cache_path <- here::here("data/cache/vrg_panint_bootstrap.rds")
  backup <- file.path(args[1],"before_vrg_panint_bootstrap.rds")
  if (file.exists(backup)) stop("Backup already exists: refuse to apply twice")
  file.copy(cache_path,backup)
  sizes <- read.csv(here::here("data/ALL_sizefreq_2024.csv")) %>%
    filter(classcode=="PANINT",campus %in% c("UCSB","VRG"),!site %in% EXCLUDED_REFERENCE_SITES)
  sites <- read.csv(here::here("data/PISCO/master_site_table_Emilyedit.csv"))
  sites <- sites[!duplicated(sites[c("site","CA_MPA_Name_Short")]),]
  pools <- merge(sizes,sites,by=c("site","site_new","campus"),all.x=TRUE) %>%
    group_by(CA_MPA_Name_Short,year,site_status) %>% group_modify(~ {
      x <- rep(.x$size,.x$count);x <- as.numeric(na.omit(x))
      data.frame(pool_n=length(x),singleton_mass=if(length(x)==1L) bio_lobster(x) else NA_real_)
    }) %>% ungroup()
  cached <- readRDS(cache_path)
  candidates <- cached %>% mutate(cache_row=row_number()) %>%
    left_join(sites %>% select(site,CA_MPA_Name_Short,site_status),by="site",relationship="many-to-many") %>%
    left_join(pools,by=c("CA_MPA_Name_Short","year","site_status")) %>%
    filter(pool_n==1L,count>0,is.finite(singleton_mass)) %>%
    transmute(cache_row,site,year,zone,count,pool_n,before_biomass=biomass,
      after_biomass=count*singleton_mass) %>% distinct()
  stopifnot(!anyDuplicated(candidates$cache_row),nrow(candidates)>0)
  cached$biomass[candidates$cache_row] <- candidates$after_biomass
  saveRDS(cached,cache_path)
  write.csv(candidates,file.path(args[1],"singleton_cache_corrections.csv"),row.names=FALSE)
  print(candidates)
})
