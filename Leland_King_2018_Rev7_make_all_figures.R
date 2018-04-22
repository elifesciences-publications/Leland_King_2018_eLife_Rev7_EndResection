

# README ------------------------------------------------------------------

# this script generates the figures from Leland ... King eLife 2018 from
# the raw data. 

# TODO: overview of directories/input data

# This code was written and tested under the following conditions:

# > sessionInfo()
# R version 3.4.4 (2018-03-15)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS High Sierra 10.13.3
# 
# Matrix products: default
# BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] grid      stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] bindrcpp_0.2.2  dplyr_0.7.4     purrr_0.2.4     readr_1.1.1     tidyr_0.8.0     tibble_1.4.2    ggplot2_2.2.1   tidyverse_1.2.1
# [9] forcats_0.3.0   lubridate_1.7.4 stringr_1.3.0   ggthemes_3.4.2  broom_0.4.4     readxl_1.1.0   
# 
# loaded via a namespace (and not attached):
#   [1] Rcpp_0.12.16        lattice_0.20-35     utf8_1.1.3          assertthat_0.2.0    digest_0.6.15       psych_1.8.3.3      
# [7] R6_2.2.2            cellranger_1.1.0    plyr_1.8.4          backports_1.1.2     acepack_1.4.1       httr_1.3.1         
# [13] pillar_1.2.1        rlang_0.2.0         lazyeval_0.2.1      rstudioapi_0.7      data.table_1.10.4-3 rpart_4.1-13       
# [19] Matrix_1.2-12       checkmate_1.8.5     splines_3.4.4       foreign_0.8-69      htmlwidgets_1.2     munsell_0.4.3      
# [25] compiler_3.4.4      modelr_0.1.1        pkgconfig_2.0.1     base64enc_0.1-3     mnormt_1.5-5        htmltools_0.3.6    
# [31] tidyselect_0.2.4    nnet_7.3-12         gridExtra_2.3       htmlTable_1.11.2    Hmisc_4.1-1         crayon_1.3.4       
# [37] nlme_3.1-131.1      jsonlite_1.5        gtable_0.2.0        magrittr_1.5        scales_0.5.0        cli_1.0.0          
# [43] stringi_1.1.7       reshape2_1.4.3      latticeExtra_0.6-28 xml2_1.2.0          Formula_1.2-2       RColorBrewer_1.1-2 
# [49] tools_3.4.4         glue_1.2.0          hms_0.4.2           parallel_3.4.4      survival_2.41-3     yaml_2.1.18        
# [55] colorspace_1.3-2    cluster_2.0.6       rvest_0.3.2         knitr_1.20          bindr_0.1.1         haven_1.1.1 



# Packages ----------------------------------------------------------------
loadNamespace("Hmisc") # for wtd.mean and wtd.var only
library(readxl)
library(broom)
library(ggthemes)
#library(ggbeeswarm)
library(stringr)
library(lubridate)
library(forcats)
library(grid)
library(tidyverse)

# ----



# ******************************************************************* -----
# HELPER FUNCTIONS, DEFAULTS, AND PLOTTING THEMES! ------------------------
# ******************************************************************* -----



# Functions and defaults related to plotting ------------------------------

# SAVE all of the figures here:
fig_save_dir <- "OUTPUT_final_pdf_versions_of_all_plots/"

# crossbar function for making semi-transparent red bars for long-range medians
crossbar_fun <-  function(x,bar_width) {
  ans <- median(x)
  tibble(y = ans, ymin = ans-(bar_width/2), ymax = ans+(bar_width/2))
}

# ggplot defaults:
gg <- list()

# vector used to add delta symbols
gg$genotype.labels <- c("wt" = expression(paste("WT")),
                        "exo1" = expression(paste(italic("exo1"),Delta)),
                        "rqh1" = expression(paste(italic("rqh1"),Delta)),
                        "crb2" = expression(paste(italic("crb2"),Delta)),
                        "crb2exo1" = expression(paste(italic("crb2"),Delta,italic("exo1"),Delta)),
                        "crb2rqh1" = expression(paste(italic("crb2"),Delta,italic("rqh1"),Delta)),
                        "rev7" = expression(paste(italic("rev7"),Delta)),
                        "rev7exo1" = expression(paste(italic("rev7"),Delta,italic("exo1"),Delta)),
                        "rev7rqh1" = expression(paste(italic("rev7"),Delta,italic("rqh1"),Delta)),
                        "rev7crb2" = expression(paste(italic("rev7"),Delta,italic("crb2"),Delta)),
                        "pht1" = expression(paste(italic("pht1"),Delta)),
                        "rev3" = expression(paste(italic("rev3"),Delta)),
                        "crb22AQ" = expression("crb22AQ"),
                        "ctp1" = expression(paste(italic("ctp1"),Delta)))

# general theme settings:
gg$my.theme <- theme_classic(22) + 
  theme(line = element_line(size=1.5, lineend = "square", color = "black"),
        axis.line = element_line(size=1.5, lineend = "square", color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        strip.background = element_blank(),
        plot.title = element_text(size = rel(1.1)))
gg$rotate.x.lab <- theme(axis.text.x = element_text(size = rel(1.2), angle = 45, hjust = 1, color = "black"))

# color settings:
# based on scale_color_tableau("tableau20")  # show_col(tableau_color_pal('tableau20')(20))
gg$colors.genotypes <- c("wt"="black",
                         "exo1"="#D62728", "rqh1"="#FF7F0E", # red and orange
                         "crb2"="#2CA02C","crb2exo1"="#98DF8A","crb2rqh1"="#BCBD22",  # greens
                         "rev7"="#9467BD","rev7exo1"="#F7B6D2","rev7rqh1"="#C5b0D5",  # purples
                         "rev7crb2"="#AeC7E8", # light blue
                         "rev3"="#7F7F7F")     # grey

# Statistics Helper Functions

# from http://stats.stackexchange.com/questions/25895/computing-standard-error-in-weighted-mean-estimation
#  Computes the variance of a weighted mean following Cochran 1977 definition
weighted_var_se <- function(x, w, na.rm=FALSE) {
  if (na.rm) { w <- w[i <- !is.na(x)]; x <- x[i] }
  # added by me: I think that sum(w) must be equal to 1:
  if (sum(w) != 1) {
    w <- w/sum(w)
  }
  n = length(w)
  xWbar = weighted.mean(x,w,na.rm=na.rm)
  wbar = mean(w)
  out = n/((n-1)*sum(w)^2)*(sum((w*x-wbar*xWbar)^2)-2*xWbar*sum((w-wbar)*(w*x-wbar*xWbar))+xWbar^2*sum((w-wbar)^2))
  return(out)
}

# Functions and defaults related to calculating resection rates -----------

# Distance between HOcs and the END of the LacO array:
kb.to.LacO.end <- 3.57+10.3

# Long-Range Resection Rate Functions

# Exclude cells and calculate rates for long-range resection rate plots:
#   - Exclude cells whose parents have a DSB (for now anyway??)
#   - Exclude cells that are sick or dead
#   - Exclude cells that do not have on-target Rad52 loading (Rad52 at frame 1 cells ARE included!)
#   - Exclude cells that load Rad52, but it happens after f.last_good_frame, or at an unknown frame (c.rad52_on = NA)
#   - calculate res.start frame ... will be 1 if Rad52 is already there at start of movie (c.rad52_on = -Inf)
#   - calculate res.end frame ... will be 1 if lacO is already gone at start of movie
#                             ... will be f.last_good_frame if it is after f.last_good_frame
#   - also make res.exact, a column that will only be TRUE if BOTH the start and end of res are exactly determined.
find_res_rates <- function (df) {
  df.res.rates <- df %>% filter(f.is_usable, 
                                is.na(c.notes) | !str_detect(c.notes, "sick|dead"),
                                !c.parent_has_dsb,
                                !is.na(c.rad52_on),
                                c.rad52_on < f.last_good_frame) %>% 
    mutate( res.start = ifelse(c.rad52_on == -Inf, 1, c.rad52_on),
            res.end  = case_when(is.na(.$c.lacO_off)                 ~ NA_real_,
                                 .$c.lacO_off == -Inf                ~ 1,
                                 .$c.lacO_off > .$f.last_good_frame  ~ .$f.last_good_frame,
                                 TRUE                                ~ .$c.lacO_off),
            res.exact = is.finite(c.rad52_on) & is.finite(c.lacO_off) & c.lacO_off <= f.last_good_frame,
            res.rate = kb.to.LacO.end / ( (res.end-res.start) * (m.t_int/60/60) )  )
  return(df.res.rates)
}

# Functions and defaults related to qPCR ----------------------------------

# Read in .xlsx files with all raw data from one or more Bio-Rad qPCR runs
# and make a list of just a few tidy data frames: run.info, summary, cycle, and melt
qpcr_import <- function(...) {
  
  # validate an even number of inputs:
  inp_args <- list(...)
  if ( length(inp_args) == 0 || 
       (length(inp_args) %% 4) != 0 ) {
    stop("Number of inputs must be even. Give inputs in groups of four: date1,plate1,xlsx.dir1,well.info.file1, date2,plate2,xlsx.dir2,well.info.file2, date3,...")
  }
  # convert to more readable table:
  inp_args <- inp_args %>% unlist() %>% matrix(nrow=4,dimnames=list(c("date","plate","xlsx.dir","well.info.file"))) %>% t() %>% as_tibble()
  
  # define this large helper function to parse a single qPCR plate
  qpcr_import_one_plate <- function(xlsx.dir,well.info.file) {
    
    # files in the xlsx.dir that will be analyzed
    c.file.contains <- c(quant.sum="Quantification Summary",end.points="End Point Results",melt.peaks="Melt Curve Peak Results",
                         quant.amp="Quantification Amplification Results",melt.amp="Melt Curve Amplification Results",melt.der="Melt Curve Derivative Results")
    c.file.endswith <- ".xlsx"
    c.file.startswith <- "^\\w.*"
    
    # excepcted columns to find in the user-generated well.info.file .xlsx file:
    c.expected.well.info <- c("well","strain","DNA","time","conc","digest","is.std","primers","rep")
    
    # initialize the FINAL output, a list called "qp"
    qp <- list()
    
    
    # find and validate all the correct .xlsx qPCR data files
    files.all  <- dir(path = xlsx.dir)
    # verifications (must match ALL patterns exactly 1 time)
    idx.matrix <- outer(files.all, str_c(c.file.startswith,c.file.contains,c.file.endswith),str_count)
    if (  any(colSums(idx.matrix) != 1) | any(rowSums(idx.matrix) > 1)  ) {
      stop(sprintf("Could not match ALL apporpriate files exactly once in the folder: \n\"%s\"",xlsx.dir))
    }
    # make a char vector of all the file names in the same order as "c.file.list"
    ind.files.keep <- apply(idx.matrix, 2, function(x)which(x==1))
    files.qpcr <- files.all[ind.files.keep]
    names(files.qpcr) <- names(c.file.contains) # copy over the same names
    
    
    
    # start by loading in the user-generated well info excel sheet
    well.info <- read_excel(file.path(xlsx.dir,well.info.file),sheet = 1)
    
    # verify that ONLY the specified col names are present!
    if (  length(setdiff(c.expected.well.info,names(well.info)) != 0)  ) {
      stop(sprintf("The column names in \"%s\" must be:\n\"%s\"",well.info.file,str_c(c.expected.well.info,collapse="\",\"")))
    }
    
    
    # helper function to convert well column to always look like this: A01,A01,B01,B02,...
    letter_numb_pad_zeros <- function(x) {
      # split the letter and the numbers into two columns
      x.split <- str_extract_all(x,"[^\\d\\s]+|\\d+",simplify=TRUE)
      x.split <- as.data.frame(x.split, stringsAsFactors = FALSE)
      
      # second column has number. Pad with zeros:
      x.split[,2] <- as.numeric(x.split[,2])      # convert to numeric
      x.dig <- floor(log10(max(x.split[,2]))) + 1 # find digits in max value
      x.split[,2] <- str_pad(x.split[,2],width=x.dig,pad="0")     # use str_pad to pad with zeros basded on digits
      
      # join the letter and numbers back together
      x <- str_c(x.split[,1],x.split[,2])
      return(x)
    }
    
    
    # clean up well.info by making chr-->factor, et.x
    well.info <- 
      well.info %>%
      mutate(well    = as.factor(letter_numb_pad_zeros(well)),
             strain   = as.factor(strain),
             digest   = as.factor(digest),
             DNA      = as.factor(DNA),
             is.std   = as.logical(is.std),
             primers  = as.factor(primers),
             rep      = as.integer(rep),
             is.noDNA = conc %in% 0,
             is.empty = is.na(strain) & is.na(DNA) & is.na(time) & is.na(conc) & is.na(digest) & is.na(is.std) & is.na(primers) & is.na(rep))  
    # now, verify that there are NO NA values anywhere in the data.frame
    # EXCEPT for in the rows that are empty wells or have no DNA
    if ( any(is.na(  well.info %>% filter(!is.noDNA,!is.empty)  )) ) {
      stop("There cannot be any NA (empty) values in the user-generated well info .xlsx except for totally empty wells")
    }
    
    
    # get summary of run information from "quant.sum" (in the 2nd excel sheet!!)
    qp$run.info <- read_excel(file.path(xlsx.dir,files.qpcr["quant.sum"]), sheet=2, col_names=FALSE)
    
    
    # build up the qp$summary df, which will have only one row per well, so it will contain
    # all data that is one per well (like Cq, end point, Melt Peak, etc.)
    tmp.cq <- read_excel(file.path(xlsx.dir,files.qpcr["quant.sum"]),sheet=1) # get well, fluor and Cq values from "Quantification Summary"
    tmp.cq <- tmp.cq %>% select(well=Well,fluor=Fluor,Cq) %>% 
      mutate(well = as.factor(letter_numb_pad_zeros(well)),
             fluor = as.factor(fluor))
    
    tmp.endpts <- read_excel(file.path(xlsx.dir,files.qpcr["end.points"]),sheet=1) # get end.RFU from "End Point Results"
    tmp.endpts <- tmp.endpts %>% select(well=Well,end.RFU=`End RFU`) %>% 
      mutate(well = as.factor(letter_numb_pad_zeros(well)))
    
    tmp.melt <- read_excel(file.path(xlsx.dir,files.qpcr["melt.peaks"]),sheet=1) # get melt-related info from "Melt Curve Peak Results"
    tmp.melt <- tmp.melt %>% select(well=Well, melt.temp=`Melt Temperature`, melt.height=`Peak Height`, melt.begin=`Begin Temperature`, melt.end=`End Temperature`) %>% 
      mutate(well = as.factor(letter_numb_pad_zeros(well)))
    
    # combine all of these together here:
    qp$summary <- 
      left_join(well.info,tmp.cq, by="well") %>% # merge with user-generated well.info
      left_join(tmp.endpts, by="well") %>%       # add in the end point data
      left_join(tmp.melt, by="well") %>%         # add in the melt info
      mutate_at(vars(melt.temp, melt.height, melt.begin, melt.end),funs(as.numeric)) %>%  # these columns have "none" in them, incorrectly interpreted as chr. Convert to dbl
      mutate(well = as.factor(letter_numb_pad_zeros(well))) # make sure well is a factor and has padded zeros
    
    
    # Make qp$cycle (has FULL traces for each amplification curve ... many rows of data per well!)
    tmp.cycle <- read_excel(file.path(xlsx.dir,files.qpcr["quant.amp"]),sheet=1)
    qp$cycle <- 
      tmp.cycle %>%
      gather(well,RFU,-Cycle) %>%
      mutate(well = as.factor(letter_numb_pad_zeros(well))) %>%
      select(well, RFU, cycle = Cycle) %>%
      left_join(qp$summary,by="well") %>%
      select(everything(),-cycle,-RFU,cycle,RFU)
    
    # Make qp$melt (has FULL traces for each melt curve ... many rows of data per well!)
    tmp.meltamp <- read_excel(file.path(xlsx.dir,files.qpcr["melt.amp"]),sheet=1) # ger raw melt curve data (RFU value at each temp for each well)
    tmp.meltamp <- tmp.meltamp %>%
      gather(well,RFU,-Temperature) %>%
      mutate(well = as.factor(letter_numb_pad_zeros(well))) %>%
      select(well, temp = Temperature, RFU)
    
    tmp.meltder <- read_excel(file.path(xlsx.dir,files.qpcr["melt.der"]),sheet=1) # get derrivative of melt curve
    tmp.meltder <- tmp.meltder %>%
      gather(well,dRFU,-Temperature) %>%
      mutate(well = as.factor(letter_numb_pad_zeros(well))) %>%
      select(well, temp = Temperature, dRFU)
    
    qp$melt <- 
      left_join(tmp.meltamp,tmp.meltder,by=c("well","temp")) %>%
      left_join(qp$summary,by="well") %>%
      select(everything(),-temp,-RFU,-dRFU,temp,RFU,dRFU)
    
    # finally, we can return qp as this fcn output!
    return(qp)
    
  }
  
  
  # loop through the rest of the rows from r=2 to the end of inp_args:
  qp <- list(run.info = tibble(), summary = tibble(), cycle = tibble(), melt = tibble())# initialize empty tibbles
  for (r in 1:nrow(inp_args)) {
    this.qp <- qpcr_import_one_plate(inp_args$xlsx.dir[r], inp_args$well.info.file[r])
    
    # add a date and plate columns to left of all data frames:
    this.qp$run.info <- this.qp$run.info %>% mutate(date=inp_args$date[r], plate=inp_args$plate[r]) %>% select(date,plate,field=1,value=2)
    this.qp$summary  <- this.qp$summary  %>% mutate(date=inp_args$date[r], plate=inp_args$plate[r]) %>% select(date,plate,everything())
    this.qp$cycle    <- this.qp$cycle    %>% mutate(date=inp_args$date[r], plate=inp_args$plate[r]) %>% select(date,plate,everything())
    this.qp$melt     <- this.qp$melt     %>% mutate(date=inp_args$date[r], plate=inp_args$plate[r]) %>% select(date,plate,everything())
    
    # now merge this.qp with the total qp (combines all plates)
    qp$run.info <- bind_rows(qp$run.info,this.qp$run.info)
    qp$summary  <- bind_rows(qp$summary,this.qp$summary)
    qp$cycle    <- bind_rows(qp$cycle,this.qp$cycle)
    qp$melt     <- bind_rows(qp$melt,this.qp$melt)
  }
  
  
  # return anything (ex: primers or digest) that has been coerced to chr from factor back to factor again:
  # NOTE: have to use this if statements bc of bug in dplyr v 0.5.0 (currently fixed in 0.5.0.9000)
  # see https://github.com/hadley/dplyr/issues/1989
  if (any(sapply(qp$summary,is.character))) {
    qp$summary <- qp$summary %>% mutate_if(is.character,as.factor)
  }
  if (any(sapply(qp$cycle,is.character))) {
    qp$cycle <- qp$cycle %>% mutate_if(is.character,as.factor)
  }
  if (any(sapply(qp$melt,is.character))) {
    qp$melt <- qp$melt %>% mutate_if(is.character,as.factor)
  }
  
  return(qp)
  
}

# Functions and defaults related to % cells w/ Rad52 foci (Fig1s3C & J) --------

# Pct Rad52-related Functions

# Exclude cells for pct Rad52-related plots:
#   - Exclude: sick, dead, out of frame, parent has on-target DSB, nuclear birth is before last good frame
#   - Exclude: the "p" cells which did nuc div just before the movie started and fission right away (c.div_nuc == -Inf)
#   - Exclude: nuclear birth timing is not excact (NA) ... this probably means it is after the last good frame anyway
#   - INCLUDE all cells where nuc birth happened before the movie started (Rad52 at T=0 counts as a DSB)
#     + this is less conservative, but the noise is huge if I don't include these cells bc my movies 
#       aren't long enough and don't start early enough to capture complete cell cycles
#   - There must be at least `min.time` minutes of usable movie after the nuclear birth (or start of movie)
#   - Do NOT include any cells that make rad52 focus, but don't meet the above critera
find_pct_rad52_valid_cells <- function (df, min.time) {
  if (missing(min.time)) stop("Please specify min.time")
  df.rad52.valid <- df %>% filter(f.is_usable,
                                  is.na(c.notes) | !str_detect(c.notes, "sick|dead"),
                                  is.na(c.notes) | !str_detect(c.notes, "out of frame"),
                                  !c.parent_has_dsb,
                                  !is.na(c.born_nuc),  # if born_nuc is NA, then we don't know if it is before or after f.last_good_frame (prob after), so we can't use it!
                                  c.born_nuc < f.last_good_frame, # c.born_nuc must be before f.last_good_frame, otherwise we def. can't use this cell
                                  !(c.div_nuc %in% -Inf)) %>% # cannot have c.div_nuc be BEFORE the movie starts (these are the dummy "p's" where fission happens right after movie starts)
    # figure out the number of minutes we've observed cell cycles that run up against the beginning
    # or end of the movie (only partially observed)
    mutate(c.is_cycle_complete = is.finite(c.born_nuc) & is.finite(c.div_nuc) & (c.div_nuc <= f.last_good_frame),
           # nuclear birth is 1 if it happens before the movie:
           c.partial_cycle_begins   = if_else(c.born_nuc == -Inf, 1, c.born_nuc), 
           # nuc division is f.last_good_frame if it happens too late or is NA (which means it prob happend after f.last_good_frame):
           c.partial_cycle_ends     = if_else(is.na(c.div_nuc) | c.div_nuc >= f.last_good_frame, f.last_good_frame, c.div_nuc)) %>%
    filter( c.is_cycle_complete |  
              ( min.time <= ((c.partial_cycle_ends-c.partial_cycle_begins)*m.t_int/60) )   
    )
  return(df.rad52.valid)
}

# Calculate Pct Rad52 from df.rad52.valid (by run and by genotype)
# TODO: need to double check how to do all these error bars!!!
calc_rad52_pcts <- function(df.rad52.valid,min.cells.cutoff) {
  if (missing(min.cells.cutoff)) stop("Please specify min.cells.cutoff")
  pct.rad52.by.run <- df.rad52.valid %>%
    # calculate % Rad52 on a per-experiment (biological replicate) basis
    group_by(e.genotype,e.folder) %>% 
    summarize(e.start.min = first(e.t_start)/60, 
              e.end.min = (first(e.t_start) + first(m.t_num)*first(m.t_int))/60, 
              e.n_valid_cells = n(), 
              e.n_Rad52_cells = sum(is.na(c.rad52_on) | is.finite(c.rad52_on) | c.rad52_on == -Inf), # count Rad52 at frame t=0!!!
              e.pct_Rad52 = e.n_Rad52_cells/e.n_valid_cells) %>%
    arrange(e.genotype,desc(e.n_valid_cells)) %>% 
    ungroup() %>% 
    # now calculate % Rad52 on a per-GENOTYPE basis
    group_by(e.genotype) %>%
    mutate(g.n_reps = n(),
           g.n_pooled_cells = sum(e.n_valid_cells),
           g.n_pooled_Rad52 = sum(e.n_Rad52_cells),
           # simple mean: just take % from each replicate and average them
           g.simple_mean = mean(e.pct_Rad52),
           g.simple_sem = sd(e.pct_Rad52)/sqrt(first(g.n_reps)),
           g.simple_ci95 = qnorm(.975) * g.simple_sem,
           # simple mean, but with cutoff >= min.cells.cutoff: just take % from each replicate and average them
           g.n_cutoff_reps  = sum(e.n_valid_cells >= min.cells.cutoff),
           g.n_cutoff_cells = sum(  ifelse(e.n_valid_cells >= min.cells.cutoff,e.pct_Rad52,NA) , na.rm=TRUE),
           g.n_cutoff_Rad52 = sum(  ifelse(e.n_valid_cells >= min.cells.cutoff,e.pct_Rad52,NA) , na.rm=TRUE),
           g.cutoff_mean    = mean( ifelse(e.n_valid_cells >= min.cells.cutoff,e.pct_Rad52,NA) , na.rm=TRUE),
           g.cutoff_sem     = sd(   ifelse(e.n_valid_cells >= min.cells.cutoff,e.pct_Rad52,NA) , na.rm=TRUE ) / sqrt(first(g.n_cutoff_reps)),
           g.cutoff_ci95 = qnorm(.975) * g.cutoff_sem,
           # weighted mean ... TODO need to check about how to do this:
           # http://stats.stackexchange.com/questions/178242/confidence-interval-for-a-weighted-mean-of-proportions
           g.wtd_mean = Hmisc::wtd.mean(e.pct_Rad52,e.n_valid_cells),
           g.wtd_sd = sqrt( Hmisc::wtd.var(e.pct_Rad52,e.n_valid_cells) ),
           g.wtd_sd_v2 = sqrt( Hmisc::wtd.var(e.pct_Rad52, e.n_valid_cells / sum(e.n_valid_cells)) ),
           g.wtd_var_se = weighted_var_se(e.pct_Rad52,e.n_valid_cells),
           g.pooled_mean = sum(e.n_Rad52_cells)/sum(e.n_valid_cells),
           g.pooled_est_se = sqrt( g.pooled_mean*(1-g.pooled_mean) / sum(e.n_valid_cells) ),
           g.pooled_ci95 = qnorm(.975) * g.pooled_est_se) %>%
    ungroup()
  return(pct.rad52.by.run)              
}

# This function counts the number of cell-cycle events in each time-bin since the
# start of HO-induction. For each bin, it finds the number of cell cycles event
# for each of these types of events (in the event.type column):
# - "nuc_begin" and "nuc_ends": # of cell cycles that being/end in that time bin, for that movie
# - "rad52_on" and "lacO_off": # of resection events that begin/end in that time bin
# - "cmpct_cyc_with_Rad52": cumulative % of cell cycles up to that time bin that have loaded Rad52
# - "cmpct_cyc_fully_resected": cumulative % of cell cycles that have fully resected lacO
# - "pct_cyc_in_resection": pct of cumulative cell cycles that are currently doing resection in that time bin
# - "pct_in_resection": similar, but something is not quite right about this ... TODO: figure that out?!
find_cdf_like_rad52_cell_counts <- function(df.rad52.valid,bin.width,bin.min,bin.max) {
  
  # define the bin edges
  bins.breaks <- seq(bin.min,bin.max,bin.width)
  # save bin info as a tibble that will be merged with the final output
  bin.info <- tibble(bin.num  = seq(1,(bin.max - bin.min) / bin.width),
                     bin.left  = seq(bin.min,bin.max-bin.width,bin.width),
                     bin.right = seq(bin.min+bin.width,bin.max,bin.width))
  
  # Note: right now, if bleaching happens, it's possible for c.partial_cycle_ends to actually be EARLIER than
  # c.rad52_on ... not sure if this matters for my analysis ... it only happens in ~20 cases.
  
  # Add new cols to df.rad52.valid that show when various things happen relative to the
  # induction of HO by adding Uracil. 
  # nuc begin/end can happen before or after the first/last frame of the movie ... if
  # that is true, then this is a "partial_cycle" but we'll just count it's begining or end
  # as happening right at e.t_start or time of the last frame of the movie
  # Also count foci that are present in the first frame of the movie as happening exactly 
  # at e.t_start (time of the first frame after Ura), but foci that happen AFTER the end
  # of the movie don't count at all ... they will still be +Inf for c.rad52_on_sec or c.lacO_off_sec
  counts.and.cumulative.counts.by.bin <- 
    df.rad52.valid %>% 
    # remove a lot of rows that we won't need (they would be dropped later on anyway)
    select(c.xlsx_row,f.order_dec16,e.genotype,e.folder,e.t_start,m.t_int,
           c.numb,c.line,c.partial_cycle_begins,c.partial_cycle_ends,c.rad52_on,c.lacO_off) %>%
    # find the time after Ura addition that these events happen for each row/cell (in minutes)
    # allow before/after movie events to count for nuc begin/end
    # only allow BEFORE movie events to count for foci events ... after movie events are 
    # left as +Inf and will end up getting dropped (become NA) during binning with cut()
    mutate(nuc_begin = ( (e.t_start + (c.partial_cycle_begins-1) * m.t_int) / 60),
           nuc_ends  = ( (e.t_start + (c.partial_cycle_ends-1)   * m.t_int) / 60),
           rad52_on  = if_else( c.rad52_on == -Inf, nuc_begin, (e.t_start + ((c.rad52_on-1)  * m.t_int)) / 60 ),
           lacO_off  = if_else( c.lacO_off == -Inf, nuc_begin, (e.t_start + ((c.lacO_off-1)  * m.t_int)) / 60 )) %>%
    # switch to "long" format: gather these four new columns into just two columns ("event.type","event.time")
    gather(key="event.type",value="event.time",nuc_begin,nuc_ends,rad52_on,lacO_off,factor_key=TRUE) %>% 
    arrange(c.xlsx_row,event.type) %>% 
    # use cut to find which bin number the event.time falls into for each row
    mutate(bin.num = cut(event.time, breaks=bins.breaks, right=TRUE, include.lowest=TRUE, labels=FALSE)) %>%
    # use count() to find out how many cells fall into each bin, then use
    # complete() to explicitly have a row for every bin, even bins with n=0
    # the number of rows is now equal to (# of exp runs ... e.folder) x (# of bin.breaks-1) x 4
    group_by(e.genotype,e.folder,event.type) %>%
    count(bin.num) %>% rename(count = n) %>%
    complete(bin.num = 1:(length(bins.breaks)-1), fill = list(count=0)) %>%
    # get the cumulative counts for each exp run (e.folder):
    mutate(e.cm.value = cumsum(count)) %>%
    ungroup() %>%
    # get the cumulative counts for each genotype
    group_by(e.genotype,event.type,bin.num) %>%
    mutate(g.cm.value = sum(e.cm.value)) %>%
    ungroup()
  
  # Somewhat tricky-er to get event.type == pct_with_Rad52 for e.cm.value and g.cm.value
  # I couldn't figure out a way to do this in a super traditionaly dplyr-y way ...
  pct.rad52.by.exp.and.genotype <- 
    counts.and.cumulative.counts.by.bin %>%
    # spread so that each event.type has its own column for e.cm.values ... have
    # to drop other columns first
    select(-count,-g.cm.value) %>%
    spread(event.type,e.cm.value) %>%
    # We can directly calculate pct_with_Rad52 (for each exp, e.__) and other related
    # pcts by exp without any grouping:
    mutate(cmpct_cyc_with_Rad52     = rad52_on/nuc_begin,
           cmpct_cyc_fully_resected = lacO_off/nuc_begin,
           pct_cyc_in_resection     = (rad52_on-lacO_off)/nuc_begin,
           pct_in_resection         = (rad52_on-lacO_off)/(nuc_begin-nuc_ends)) %>%
    # for values that are cumulative per genotype (combining all exps), we need
    # to group first, then use sum() (note: don't need cumsum, because these
    # values are already cumulative sums on a per exp basis)
    group_by(e.genotype,bin.num) %>%
    mutate(g.cm.cmpct_cyc_with_Rad52     = sum(rad52_on) / sum(nuc_begin),
           g.cm.cmpct_cyc_fully_resected = sum(lacO_off) / sum(nuc_begin),
           g.cm.pct_cyc_in_resection     = (sum(rad52_on) - sum(lacO_off)) / sum(nuc_begin),
           g.cm.pct_in_resection         = (sum(rad52_on) - sum(lacO_off)) / (sum(nuc_begin) - sum(nuc_ends))) %>%
    ungroup() %>%
    # now gather back the e.___ values into two columns, and select out only the cols 
    # that will be appended onto rad52.cell.counts that we've already made above.
    gather(key = "event.type", value = "e.cm.value", nuc_begin:pct_in_resection, factor_key = TRUE) %>%
    select(e.genotype,e.folder,bin.num,event.type,e.cm.value,everything()) %>%
    # NOTE: at this point we INCORRECTLY have a g.cm.___ values that repeats the 
    # cumulative over genotype values for over and over, even when
    # event.type != pct_... or cmpct_... So, drop all rows except for the newly
    # created pct ones
    filter(event.type %in% c("cmpct_cyc_with_Rad52","cmpct_cyc_fully_resected","pct_cyc_in_resection","pct_in_resection")) %>%
    # and finally, this is also kind of hack-y, but use a major cases_when statement to combine all of the
    # g.cm.___ columns into a single column with the appropriate values for each row based on
    # the event.type for that row:
    mutate(g.cm.value = case_when(
      .$event.type == "cmpct_cyc_with_Rad52"     ~ .$g.cm.cmpct_cyc_with_Rad52,
      .$event.type == "cmpct_cyc_fully_resected" ~ .$g.cm.cmpct_cyc_fully_resected,
      .$event.type == "pct_cyc_in_resection"     ~ .$g.cm.pct_cyc_in_resection,
      .$event.type == "pct_in_resection"         ~ .$g.cm.pct_in_resection
    )) %>%
    select(-(g.cm.cmpct_cyc_with_Rad52:g.cm.pct_in_resection))
  
  
  # now use bind_rows to combine these two together, giving us all event.types
  # and their counts (if applicable) and cumulative sums/pcts by exp (e.cm.value) 
  # and cumulative sums/pcts by genotype (g.cm.value)
  rad52.cumulative.counts <- 
    bind_rows(counts.and.cumulative.counts.by.bin, pct.rad52.by.exp.and.genotype) %>% 
    mutate(event.type = factor(event.type,levels=c("nuc_begin","nuc_ends","rad52_on","lacO_off",
                                                   "cmpct_cyc_with_Rad52","cmpct_cyc_fully_resected",
                                                   "pct_cyc_in_resection","pct_in_resection"))) %>%
    left_join(bin.info) %>%
    select(e.genotype,e.folder,event.type,bin.num,bin.left,bin.right,count,e.cm.value,g.cm.value) %>%
    arrange(e.genotype,e.folder,event.type,bin.num)
  
  # this is the function's output:
  return(rad52.cumulative.counts)
}


# ----



# ******************************************************************* -----
# IMPORT THE DATA HERE! ---------------------------------------------------
# ******************************************************************* -----



# Import the sub-pixel foci intensity quantification data here ------------
# for plots in figures 1D, 1E, 1supp2A', 1supp2B', 2B, etc.

tmp.rdata.save.file.path <- "OUTPUT_processed_data_saved_as_RData_files/all_subpixel_foci_intensity_quantification.RData"
if ( file.exists(tmp.rdata.save.file.path) ) {
  
  # load in the already processed data.frame:
  load(tmp.rdata.save.file.path)
  
} else {
  
  # generate the data.frame by combining all raw .csv files:
  
  tmp.dir <- "INPUT_raw_data/raw_subpixel_foci_intensity_data/"
  tmp.data.to.read.in <- tribble(
    ~csv_file_with_measurments,                                                              ~figure,       ~genotype, ~field, ~cell, ~t.start, ~t.end, ~t.RadOn, ~t.LacOff,
    "rev7_F2S1.1_f37c7b_Crp61_2015_0724_run1_2149rev7_60min__10_Measurments_v1.csv",         "Figure2s1A'", "rev7",    37,     "7b",  1,        31,     23,       29,        
    "rev7_F2S1.2_f335c3p_Crp61_2016_0215_run2_2149rev7_120min_14_Measurments_v1.csv",        "Figure2s1B'", "rev7",    335,    "3p",  1,        31,     3,        7,      
    "rev7_F2S1.3_f371c17b_Crp61_2015_0724_run2_2149rev7_50min__01_Measurments_v1.csv",       "Figure2s1C'", "rev7",    371,    "17b", 1,        31,     19,       23,
    "rev7_F2S1.4_f316c13p_Crp61_2016_0209_2149rev7_later_run1_80MIN_03_Measurments_v1.csv",  "Figure2s1D'", "rev7",    316,    "13p", 1,        25,     2,        6,
    "rev7_Fig2A_f320c15b_Crp61_2016_0215_run2_2149rev7_120min_12_Measurments_v1.csv",        "Figure2B",    "rev7",    320,    "15b", 1,        31,     9,        17,
    "wt_Fig1B_f61c8a_Crp61_2014_1004_run3_1771wt_50min_08_Measurments_v1.csv",               "Figure1D",    "wt",      61,     "8a",  1,        31,     7,        22,
    "wt_Fig1B_f61c8b_Crp61_2014_1004_run3_1771wt_50min_08_Measurments_v1.csv",               "Figure1E",    "wt",      61,     "8b",  1,        31,     7,        21,
    "wt_Fig1S2.1_f282c11b_Crp61_2015_0429_run2_1771_45min_19_Measurments_v1.csv",            "Figure1s2A'", "wt",      282,    "11b", 1,        31,     5,        15,
    "wt_Fig1S2.2_f98c10b_Crp61_2015_0429_run2_1771_45min_10_Measurments_v1.csv",             "Figure1s2B'", "wt",      98,     "10b", 1,        31,     4,        22
  )
  
  df_subpx <- 
    tmp.data.to.read.in %>%
    group_by_all() %>%
    do(read_csv(str_c(tmp.dir,.$csv_file_with_measurments))) %>%
    ungroup() %>%
    rename(row = X1) 
  
  # check that the proper number of rows got imported for everything (some movies are longer)
  df_subpx %>% group_by(csv_file_with_measurments) %>% tally()
  
  # save this file for faster importing later, if needed:
  save(df_subpx, file=tmp.rdata.save.file.path)
}



# Import and clean up the resection duration scoring data -----------------
# for beehive plots of resection rate, e.g. figures 1F, 2C, 3A

tmp.rdata.save.file.path <- "OUTPUT_processed_data_saved_as_RData_files/Leland2018_resection_assay_raw_scoring_data.RData"
if ( file.exists(tmp.rdata.save.file.path) ) {
  
  # load in the already processed data.frame:
  load(tmp.rdata.save.file.path)

} else {
  
  # The required file Leland2018_resection_assay_raw_scoring_data.RData is generated using a combination 
  # of the following files.
  # Eventually I may copy over that code into this section.
  
  # import_cell_data_validate.R
  # parse_dv_log_files_for_metadata.R
  # randomizing_script.R
  # flds_with_metadata.RData
  
  # Ideally, I'd write this as a helper function in the top section, to make this more readable hopefully
  
}


# Import the older qPCR data from Robert’s 2/3/16 exp (Fig1s3C) -----------

# Import the raw data from Robert's 2016_0203 qPCR experiment.
# this .RData file is basically just the raw qPCR run data converted from .xlsx 
# files to data.frames using qpcr_import() below. If this has already been
# done once, then we will just load the .RData file. If this .RData file does
# not yet exist, then we will make it in the qpcr_bargraph block of code
# by importing the .xlsx files using qpcr_import()


tmp.rdata.save.file.path <- "OUTPUT_processed_data_saved_as_RData_files/2016_0203_qPCR_WT_Robert.RData"
if ( file.exists(tmp.rdata.save.file.path) ) {
  
  # load directly from already made .RData file, if it exists
  load(tmp.rdata.save.file.path)

} else {

  # use the helper function to import the raw excel data:
  qp_old <- qpcr_import("2016_0203",1,"INPUT_raw_data/raw_qPCR_data/2016_0203_qPCR_WT_Robert/","2016_0203_well_info.xlsx")
  # save this as an .RData file to speed things up next time around:
  save(qp_old,file = tmp.rdata.save.file.path)
  
}



# Import the newer qPCR data from Megan’s exps (Fig 1s3H&I, Fig2D)  -------

# Same thing, but for all the more recent qPCR data that Megan ran:
dir_for_new_qpcr_data <- "INPUT_raw_data/raw_qPCR_data/"

tmp.rdata.save.file.path <- "OUTPUT_processed_data_saved_as_RData_files/all_newer_qPCR_data_as_of_2017_0809.RData"

if ( file.exists(tmp.rdata.save.file.path)  ) {
  
  # directly load the already made .RData file (made in else statement below) if it exists, to save time
  load(tmp.rdata.save.file.path)
  
} else {
  
  # order of inputs for each plate is: "date",plate_number,"xlsx_dir_name","well_info_file_name"
  tmp.qpcrimport.fcn.args = list(
    "2017_0628", 1, file.path(dir_for_new_qpcr_data,"2017_0628_first_try_new_PlusMinus_primers","RawData"), "2017_0628_well_info.xlsx",
    "2017_0703", 1, file.path(dir_for_new_qpcr_data,"2017_0703_megan1_plusminus13000","RawData"), "2017_0703_well_info.xlsx",
    "2017_0705", 1, file.path(dir_for_new_qpcr_data,"2017_0705_megan_DiffGenotypes_m168","RawData"), "2017_0705_well_info.xlsx",
    "2017_0706", 1, file.path(dir_for_new_qpcr_data,"2017_0706_megan_m300_p590_rev7_wt","RawData"), "2017_0706_well_info_MCK_v2.xlsx",
    "2017_0714", 1, file.path(dir_for_new_qpcr_data,"2017_0714_megan_2456rev7rqh1","RawData"), "2017_07014_well_info_BAL_fromPic.xlsx",
    "2017_0727", 1, file.path(dir_for_new_qpcr_data,"2017_0727_megan_2123wt_long_induc","RawData"), "2017_07027_well_info.xlsx",
    "2017_0728", 1, file.path(dir_for_new_qpcr_data,"2017_0728_megan_2149rev7_fill_in","RawData"), "2017_0728_well_info_MCK.xlsx",
    "2017_0807", 1, file.path(dir_for_new_qpcr_data,"2017_0807_megan_wt_rqh1_180","RawData"), "2017_0807_well_info_MCK.xlsx",
    "2017_0809", 1, file.path(dir_for_new_qpcr_data,"2017_0809_megan_exo1_short","RawData"), "2017_0809_well_info_MCK.xlsx",
    "2017_0912", 1, file.path(dir_for_new_qpcr_data,"2017_0912_megan_1500_new","RawData"), "2017_0912_well_info_MCK_BLedits.xlsx"
  )
  # use the helper function qpcr_import() from above to import the raw excel data:
  qp <- do.call(qpcr_import,tmp.qpcrimport.fcn.args)
  
  # check out data and clean up a bit if needed:
  lapply(qp,function(x) sapply(x,class))
  lapply(qp$summary, levels)
  
  # reorder ALL factor levels in qp$summary, qp$cycle, and qp$melt:
  tmp.fct.order <- list(strain = c("2123","1914","2149","2477","2183","2456"),
                        digest = c("HincII","ApoI","none"),
                        primers = c("Ncb2_v1","Ncb2_v2","m168_v2","m300","m3023","m14253","p590","p2996","p13150"))
  for (df.name in c("summary","cycle","melt")) {
    for (var.name in names(tmp.fct.order)) {
      qp[[df.name]][[var.name]] <- fct_relevel(qp[[df.name]][[var.name]], tmp.fct.order[[var.name]])
    }
  }
  
  # Use left_join to make a new variable is.expWell
  # generally these are all the non-std curve wells with 2.5 ng/uL gDNA in them, 
  # HOWEVER, sometime the 2.5 ng/uL point on the std curve is used as the HincII
  # digest experimental well, so a well can be BOTH a std curve well AND an exp well!
  tmp.is.expWells <- 
    qp$summary %>%
    filter(!is.empty, !is.noDNA) %>%
    select(date,plate,strain,time,conc,digest,primers,is.std) %>%
    distinct() %>%
    mutate(is.expWell = 
             case_when(
               .$date=="2017_0628" & .$conc==2.5 & .$digest!="none" & .$primers!="Ncb2_v2" ~ TRUE,
               .$date %in% c("2017_0703","2017_0705","2017_0706","2017_0714","2017_0727","2017_0728","2017_0807","2017_0809") & .$conc==2.5 ~ TRUE,
               .$date == "2017_0912" & .$conc == 2.5  & .$is.std==FALSE ~ TRUE,
               TRUE ~ FALSE)) %>% arrange(date,is.expWell)
  # here are the left joins, followed by rearranging the columns into a sensible order:
  for (df.name in c("summary","cycle","melt")) {
    qp[[df.name]] <- 
      left_join(qp[[df.name]], tmp.is.expWells) %>% 
      mutate(is.expWell = if_else(is.na(is.expWell), FALSE, is.expWell))
  }
  
  # Make is.excluded: a column to mark all of the excluded wells for each exp
  tmp.is.excluded <- 
    qp$summary %>%
    filter( (date == "2017_0628" & well %in% c("B04","F05","D05")) |
              (date == "2017_0706" & well == "C05") |
              (date == "2017_0705" & well %in% c("C10","C12")) |
              (date == "2017_0714" & well %in% c("C11","D11")) |
              (date == "2017_0728" & well %in% c("C02","B03")) |
              (date == "2017_0809" & well %in% c("F11","B10")) |
              (date == "2017_0912" & well %in% c("A04","B04"))
    ) %>%
    select(date,plate,well) %>%
    mutate(is.excluded = TRUE)
  
  # here are the left joins, followed by rearranging the columns into a sensible order:
  for (df.name in c("summary","cycle","melt")) {
    qp[[df.name]] <- 
      left_join(qp[[df.name]], tmp.is.excluded) %>% 
      mutate(is.excluded = if_else(is.na(is.excluded), FALSE, is.excluded))
  }
  
  # reorder columns appropriately for qp$summary, qp$cycle, and qp$melt:
  for (df.name in c("summary","cycle","melt")) {
    qp[[df.name]] <- 
      qp[[df.name]] %>%
      select(date,plate,well,strain,DNA,time,conc,digest,primers,rep,is.empty,is.noDNA,is.std,is.expWell,is.excluded,everything())
  }
  
  # save this as an .RData file to speed things up next time around:
  save(qp, file = tmp.rdata.save.file.path )
}


# now that this `pq` data.frame is imported, perform some additiona
# calculations to make the std curves and linear fits that will allow
# conversion from raw CQ values to DNA concentrations: (for subsequent
# plotting)
qp$lin.fits <- 
  qp$summary %>% 
  filter(is.std, !is.excluded) %>% 
  group_by(date,plate,primers,digest,strain,time) %>% # these are all the col variables that distinguish the std curves.
  do(fit = lm(Cq ~ log(conc), data = .) ) %>%  # fit a line to each std curve using lm()
  # now this gets ugly, but basically this is pulling out the y.intercept, slope, and R^2 values
  # by using both the tidy() and glance() functions from broom, and then combining them 
  # together with a left_join
  { left_join(
    tidy(.,fit) %>%   # tidy pulls out the slope and y.intercept, but we have to use spread() to make them separate columns
      select(primers,digest,term,estimate) %>% 
      spread(term,estimate) %>% rename(y.intercept=`(Intercept)`, slope=`log(conc)`),
    glance(.,fit) %>% # glance pulls out the r.sqrd values 
      select(primers,digest,r.squared, adj.r.squared)
  )
  } %>% 
  ungroup() %>%
  mutate(  qpcr.effic = exp(-1/slope) - 1  ) # qPCR efficiency = e^(-1/slope) - 1

# view the linear fits:
qp$lin.fits

# manually select out ONLY the unique linear fits (one for each primer) that I will use:
qp$lin.fits <- 
  qp$lin.fits %>%
  mutate(used.for.fits = digest != "none" & !(primers == "Ncb2_v1" & date == "2017_0703"))

# Use linear fits to calculate concentration (calc.conc) from Cq:
qp$summ.calc.conc <- 
  qp$lin.fits %>%
  # only merge in unique linear fits manually selected above
  filter(used.for.fits) %>% 
  arrange(primers) %>%
  # do left_join based ONLY on primers
  select(primers,y.intercept:length(.)) %>%
  left_join(qp$summary,.) %>% 
  mutate(calc.conc = case_when(.$is.empty  ~ NA_real_, # empty wells have a calc.conc of NA
                               is.na(.$Cq) ~ 0,        # all other wells with Cq=NA have a calc.conc of 0
                               TRUE        ~ exp( (.$Cq-.$y.intercept)/.$slope) )) # this is the formula to convert from Cq to calc.conc w/ natural log



# Import the analyzed qPCR results from primers across HOcs (Fig 1SA, 2E) --------

hocs.qpcr <- 
  read_csv("INPUT_raw_data/raw_qPCR_data/analyzed_qPCR_results_primers_across_HOcs.csv", skip=2)


# ----



# ******************************************************************* -----
# MAKE THE PLOTS HERE! ----------------------------------------------------
# ******************************************************************* -----



# Plots of sub-pixel foci quantification data (Fig1D, 1E, 1s2A', etc) --------

# hardcode the pixel size (from image metadata)
tmp.px.size.xy.in.um <- 0.1605
tmp.px.size.z.in.um <- 0.26

# params for sub-pixel resolution (same for every plot)
tmp.interp <- "None"
tmp.inc_res_by <- 24

# these "best" params were determined emperically above
tmp.h_signal <- 2 #number of z slices above AND below to include in the "signal" cylindar
tmp.r_signal <- 3 #radius in pixels of the signal cylindar
tmp.h_bkg <- 5
tmp.r_bkg <- 4

tmp.figures <- df_subpx$figure %>% unique() # each figure to do


tmp.gg.my.theme <- theme_classic(22) + 
  theme(line = element_line(size=1.5, lineend = "square", color = "black"),
        axis.line = element_line(size=1.5, lineend = "square", color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        strip.background = element_blank(),
        plot.title = element_text(size = rel(1.1)))

# loop through and make a plot for EACH combination of these measurment params in tmp.comb.of.cutoffs
for (this_fig in tmp.figures) {
  
  this__comb.of.cutoffs <- 
    crossing(h_signal = tmp.h_signal, r_signal = tmp.r_signal,
             bkg_h_above_sig_h = tmp.h_bkg - tmp.h_signal, bkg_r_beyond_sig_r = tmp.r_bkg - tmp.r_signal) %>%
    mutate(h_bkg = h_signal + bkg_h_above_sig_h,
           r_bkg = r_signal + bkg_r_beyond_sig_r)
  
  
  tmp.toplot <- 
    df_subpx %>%
    filter(figure == this_fig,
           timepoint >= t.start & timepoint <= t.end, # crop off timepoints that are not displayed in figs
           interp == tmp.interp,  # only select the measurments that are needed for these specific conditions defined above
           inc_res_by == tmp.inc_res_by,
           abs(z_measured) <= tmp.h_bkg,
           radius_measured %in% c(tmp.r_signal,tmp.r_bkg)) %>%
    
    mutate(scale_and_interp = factor(str_c(inc_res_by,interp), 
                                     levels=c("1None","24None","24Bicubic"),
                                     labels=c("pixel resolution", "subpixel res 24x24","subpixel interp"))) %>%
    
    select(csv_file_with_measurments:t.LacOff, scale_and_interp,   # these will facet the plots later
           timepoint, channel, # these variables are for the x-axis and color
           Area, IntDen, z_measured, radius_measured) %>%  # these will be used to compute intensity for signal, bkg, and bkg subtracted
    
    # now create a ton of extra rows, specifying all the combinations of
    # signal radius, signal cylinder height, background radius, and background
    # cylinder height that you want to try. Note that cylinder "height" is the 
    # number of slices above and below the middle slice, so "height" of 1 means
    # include slices -1, 0, +1. "height of 3 means -3,-2,-1,0,1,2,3
    crossing(this__comb.of.cutoffs) %>%
    
    # now prune again. Get rid of any rows with z slices or radii outside the range:
    filter(radius_measured <= r_bkg,
           abs(z_measured) <= h_bkg) %>%
    
    group_by_at(vars(csv_file_with_measurments:t.LacOff, 
                     timepoint, channel, scale_and_interp, 
                     h_signal, r_signal, bkg_h_above_sig_h, bkg_r_beyond_sig_r, h_bkg, r_bkg)) %>%
    
    summarize(IntDen_sigCyl  = sum(IntDen[radius_measured == r_signal & abs(z_measured) <= h_signal]),
              IntDen_bkgCyl  = sum(IntDen[radius_measured == r_bkg    & abs(z_measured) <= h_bkg]),
              Area_sigCyl = first(Area[radius_measured == r_signal]) * (tmp.px.size.xy.in.um^2) *
                ((2*first(h_signal))+1) * tmp.px.size.z.in.um,
              Area_bkgCyl = first(Area[radius_measured == r_bkg]) * (tmp.px.size.xy.in.um^2) *
                ((2*first(h_bkg))+1) * tmp.px.size.z.in.um) %>%
    ungroup() %>%
    
    # now calculate average intensity and background subtracted version
    mutate(AveInt_sig = IntDen_sigCyl / Area_sigCyl,
           AveInt_bkg = (IntDen_bkgCyl - IntDen_sigCyl) / (Area_bkgCyl - Area_sigCyl),
           AveInt_sigMinusBkg = AveInt_sig - AveInt_bkg) %>%
    
    # normalize (could tweak how this is done maybe?)
    group_by_at(vars(csv_file_with_measurments:t.LacOff, 
                     channel, scale_and_interp)) %>%
    
    mutate(AveInt_normSigMinusBkg = (AveInt_sigMinusBkg - min(AveInt_sigMinusBkg)) / (max(AveInt_sigMinusBkg) - min(AveInt_sigMinusBkg)))
  
  
  
  p <- 
    tmp.toplot %>%
    ggplot(aes(y=AveInt_normSigMinusBkg, x=timepoint, color=as.factor(channel),
               group = interaction(channel,scale_and_interp))) +
    geom_vline(aes(xintercept = t.RadOn), color = "#0000FF", size = 1.8, linetype="31") + # first frame
    geom_vline(aes(xintercept = t.LacOff), color = "#0000FF", size = 1.8, linetype="31") + # last frame
    geom_line(size = 2) + 
    scale_y_continuous(name = "Normalized focus\nintensity (A.U.)") + 
    scale_x_continuous(name = expression(paste("Frame (",Delta,"t = 10 min)"))) + 
    scale_color_manual(name = "", values = c("#39A04A","#FF00FF"), labels = c("LacI-GFP","Rad52-mChr")) + 
    coord_fixed(ratio = 15) + 
    tmp.gg.my.theme + 
    theme(legend.text.align=0,
          legend.key.size = unit(2.1, 'lines'),
          legend.text = element_text(size = rel(0.8), color = "black"))
  ggsave(plot = p,
         filename = str_c(fig_save_dir,this_fig,"_subpixel_focus_intensity.pdf"),
         device = cairo_pdf(),
         width = 40, height = 9, units = "cm")
  #dev.off()
}








# Beehive plots of resection rate (Fig 1F, 2C, 3A, 3s2B) ------------------------

# calculate all resection rates from raw datausing helper function
df.res.rates <- find_res_rates(df)


# Figure 1F: beehive plot of resection rate for WT, exo1, crb2
tmp.genotypes <- c("wt","exo1","crb2")

# fancy way of breaking this up into two steps so we can get "binned" versions
# of the red median crossbars that align perfectly with the dots!
tmp.g <- 
  df.res.rates %>% filter(f.is_usable, res.exact, is.na(c.notes) | !str_detect(c.notes, "sick|dead"),
                          e.genotype %in% tmp.genotypes) %>%
  mutate(e.genotype = factor( e.genotype, levels = tmp.genotypes)) %>%
  ggplot(aes(x=e.genotype,y=res.rate, group=e.genotype)) + 
  geom_dotplot(binaxis = "y", stackdir = "center", 
               method="histodot",binwidth=0.8, dotsize=.9)
tmp.g +
  stat_summary(data=layer_data(tmp.g) %>% mutate(e.genotype=factor(x,levels=1:(length(tmp.genotypes)-1),labels=tmp.genotypes[-2])),
               aes(x = e.genotype, y = y),
               geom="crossbar", fun.data=crossbar_fun, fun.args=c(bar_width=.75),
               color=NA,fill='red',size=0.8,width=1.2,alpha=0.5) + 
  scale_x_discrete(name="", labels=gg$genotype.labels, drop=FALSE) + 
  scale_y_continuous(name="Resection rate (kb/hr)") + 
  gg$my.theme + gg$rotate.x.lab + 
  theme(panel.background = element_blank(),
        plot.background = element_blank())

ggsave(filename = str_c(fig_save_dir,"Figure1F_res_rate_wt_exo1_crb2.pdf"),
       device = cairo_pdf(),
       scale = 1, width = 11, height = 14, units = "cm")




# Figure 2C: beehive plot of resection rate for WT, exo1, crb2
tmp.genotypes <- c("wt","rev7", "crb2","rev3")

# fancy way of breaking this up into two steps so we can get "binned" versions
# of the red median crossbars that align perfectly with the dots!
tmp.g <- 
  df.res.rates %>% filter(f.is_usable, res.exact, is.na(c.notes) | !str_detect(c.notes, "sick|dead"),
                          e.genotype %in% tmp.genotypes) %>%
  mutate(e.genotype = factor( e.genotype, levels = tmp.genotypes)) %>%
  ggplot(aes(x=e.genotype,y=res.rate, group=e.genotype)) + 
  geom_dotplot(binaxis = "y", stackdir = "center", 
               method="histodot",binwidth=0.8, dotsize=.9)
tmp.g +
  stat_summary(data=layer_data(tmp.g) %>% mutate(e.genotype=factor(x,levels=1:(length(tmp.genotypes)),labels=tmp.genotypes)),
               aes(x = e.genotype, y = y),
               geom="crossbar", fun.data=crossbar_fun, fun.args=c(bar_width=.75),
               color=NA,fill='red',size=0.8,width=0.95,alpha=0.5) + 
  scale_x_discrete(name="", labels=gg$genotype.labels, drop=FALSE) + 
  scale_y_continuous(name="Resection rate (kb/hr)") + 
  gg$my.theme + gg$rotate.x.lab + 
  theme(panel.background = element_blank(),
        plot.background = element_blank())

ggsave(filename = str_c(fig_save_dir,"Figure2C_res_rate_wt_rev7_crb2_rev3.pdf"),
       device = cairo_pdf(),
       scale = 1, width = 13, height = 14, units = "cm")





# Figure 3A: beehive plot of resection rate, epistasis for all genotypes
tmp.genotypes <- c("wt","exo1","rqh1","crb2","crb2exo1","crb2rqh1","rev7","rev7exo1","rev7rqh1","rev7crb2")

# fancy way of breaking this up into two steps so we can get "binned" versions
# of the red median crossbars that match perfectly!
tmp.g <- 
  df.res.rates %>% filter(f.is_usable, res.exact, is.na(c.notes) | !str_detect(c.notes, "sick|dead"),
                          e.genotype %in% tmp.genotypes) %>%
  mutate(e.genotype = factor( e.genotype, levels = tmp.genotypes)) %>%
  ggplot(aes(x=e.genotype,y=res.rate, group=e.genotype)) + 
  geom_dotplot(binaxis = "y", stackdir = "center", 
               method="histodot",binwidth=0.8, dotsize=.9)
tmp.g +
  stat_summary(data=layer_data(tmp.g) %>% mutate(e.genotype=factor(x,levels=1:(length(tmp.genotypes)-1),labels=tmp.genotypes[-2])),
               aes(x = e.genotype, y = y),
               geom="crossbar", fun.data=crossbar_fun, fun.args=c(bar_width=.75),
               color=NA,fill='red',size=0.8,width=.8,alpha=0.5) + 
  scale_x_discrete(name="", labels=gg$genotype.labels, drop=FALSE) + 
  scale_y_continuous(name="Resection rate (kb/hr)") + 
  gg$my.theme + gg$rotate.x.lab + 
  theme(panel.background = element_blank(),
        plot.background = element_blank())

ggsave(filename = str_c(fig_save_dir,"Figure3A_res_rate_all_genotypes_epistasis.pdf"),
       device = cairo_pdf(),
       scale = 1, width = 29, height = 16, units = "cm")




# Figure 3s2B: beehive plot of resection DURATION (flipped axis) for all genotypes
# including both exactly determined durations (black) and non-exactly determined (red)

tmp.genotypes <- c("wt","exo1","rqh1","crb2","crb2exo1","crb2rqh1","rev7","rev7exo1","rev7rqh1","rev7crb2","rev3")

tmp.g <- df.res.rates %>% filter(f.is_usable, is.na(c.notes) | !str_detect(c.notes, "sick|dead"),
                                 e.genotype %in% tmp.genotypes,
                                 !(e.genotype == "exo1" & res.exact)) %>%    
  # leave an empty space to write in "exo1":
  mutate(e.genotype = factor( e.genotype, levels = tmp.genotypes ),
         res.durtn = (res.end-res.start) * (m.t_int/60)) %>%
  ggplot(aes(x=e.genotype,y=res.durtn,fill=!res.exact,color=!res.exact)) + 
  geom_dotplot(binaxis = "y", stackdir = "center", 
               method="histodot", binwidth=10, dotsize=.7, 
               position=position_dodge(width = 0.8)) +
  scale_x_discrete(name="", labels=gg$genotype.labels, drop=FALSE) +
  scale_y_continuous(name="Resection duration (min)") + 
  scale_fill_manual(values=c("black","red"), labels=c("exact duration","minimum possible duration")) + 
  scale_color_manual(values=c("black","red"), labels=c("exact duration","minimum possible duration")) +
  gg$my.theme + gg$rotate.x.lab +
  theme(panel.background = element_blank(),
        plot.background = element_blank(),
        legend.position = c(.75,.9),
        legend.direction = "vertical",
        legend.title = element_blank(),
        legend.margin = margin(0,8,6,4),
        legend.text = element_text(size = rel(0.8)),
        legend.background = element_rect(color = "black"))

ggsave(tmp.g,
       filename = str_c(fig_save_dir,"Figure3s2B_res_DURATION_all_genotypes_plus_non-exact.pdf"),
       device = cairo_pdf(),
       scale = 1, width = 31, height = 14, units = "cm")

# make the 2nd y-axis (resection rate)
decm.places <- 2
dur.to.rate <- function(x) format(round(kb.to.LacO.end/(x/60),decm.places),nsmall=decm.places)
rate.to.dur <- function(x) round((kb.to.LacO.end/x)*60,6)
tmp.color <- "gray40"
tmp.g2 <- 
  tmp.g + 
  scale_y_continuous(name="Resection rate (kb/hr)",
                     breaks=rate.to.dur(as.double(dur.to.rate(c(100,200,300,400)))),
                     labels=dur.to.rate) +
  theme(text = element_text(color=tmp.color),
        axis.line = element_line(color=tmp.color),
        axis.ticks = element_line(color=tmp.color),
        axis.text = element_text(color=tmp.color),
        rect = element_rect(color=tmp.color))

ggsave(tmp.g2,
       filename = str_c(fig_save_dir,"Figure3s2B_second_y-axis_showing_rate.pdf"),
       device = cairo_pdf(),
       scale = 1, width = 27, height = 14, units = "cm")



# qPCR with primers across the HOcs (Fig 1sA, 2E) -----------------------------------

# Fig 1sA: ONLY plot the wt and exo data
hocs.qpcr %>%
  gather(key="genotype",value="pct.cut",wt,rev7,exo1) %>%
  filter(genotype %in% c("wt","exo1")) %>%
  mutate(genotype = factor(genotype,levels = c("wt","exo1"))) %>%
  ggplot(aes(x=genotype, y=pct.cut)) +
  stat_summary(fun.y = mean, geom="bar", color="black", fill="grey") +
  stat_summary(fun.ymin = function(x) mean(x) - ((sd(x)/sqrt(length(x))) * qnorm(0.975)), 
               fun.ymax = function(x) mean(x) + ((sd(x)/sqrt(length(x))) * qnorm(0.975)), 
               geom="errorbar", width = 0.2, size = 1) +
  scale_y_continuous(name = "Percent loss of qPCR\nproduct across the HOcs", labels = function(x) str_c(x,"%")) +
  scale_x_discrete(name="", labels = gg$genotype.labels) + 
  gg$my.theme + 
  theme(plot.title = element_text(hjust=0.5))

ggsave(filename = str_c(fig_save_dir,"Figure1s3A_HOcsQPCR_WTexo1Only.pdf"),
       device = cairo_pdf(),
       scale = 1, width = 8, height = 11, units = "cm")

# stats for drawing in on this figure:
hocs.qpcr %>%
  gather(key="genotype",value="pct.cut",wt,rev7,exo1) %>%
  mutate(genotype = factor(genotype,levels = c("wt","rev7","exo1"))) %>%
  {pairwise.t.test(.$pct.cut, .$genotype, 
                   p.adjust.method = "none",
                   pool.sd = FALSE,
                   paired = FALSE,
                   alternative = "two.sided")}


# Fig 2E: ONLY plot the wt and rev7 data
hocs.qpcr %>%
  gather(key="genotype",value="pct.cut",wt,rev7,exo1) %>%
  filter(genotype %in% c("wt","rev7")) %>%
  mutate(genotype = factor(genotype,levels = c("wt","rev7"))) %>%
  ggplot(aes(x=genotype, y=pct.cut)) +
  stat_summary(fun.y = mean, geom="bar", color="black", fill="grey") +
  stat_summary(fun.ymin = function(x) mean(x) - ((sd(x)/sqrt(length(x))) * qnorm(0.975)), 
               fun.ymax = function(x) mean(x) + ((sd(x)/sqrt(length(x))) * qnorm(0.975)), 
               geom="errorbar", width = 0.2, size = 1) +
  scale_y_continuous(name = "Percent loss of qPCR\nproduct across the HOcs", labels = function(x) str_c(x,"%")) +
  scale_x_discrete(name="", labels = gg$genotype.labels) + 
  gg$my.theme + 
  theme(plot.title = element_text(hjust=0.5))

ggsave(filename = str_c(fig_save_dir,"Figure2E_HOcsQPCR_WTrev7Only.pdf"),
       device = cairo_pdf(),
       scale = 1, width = 8, height = 11, units = "cm")

# stats for drawing in on this figure:
hocs.qpcr %>%
  gather(key="genotype",value="pct.cut",wt,rev7,exo1) %>%
  mutate(genotype = factor(genotype,levels = c("wt","rev7","exo1"))) %>%
  {pairwise.t.test(.$pct.cut, .$genotype, 
                   p.adjust.method = "none",
                   pool.sd = FALSE,
                   paired = FALSE,
                   alternative = "two.sided")}









# Combined plot: qPCR @ 168nt and %cells w/ Rad52 foci (Fig 1s3C) --------

# fit lines using lm() to find y.intercept, slope, r.squared, adj.r.squared,
# and qpcr.effic for each standard curve

# for this particular exp, I fit separate std curves for each timepoint (mostly
# these curvues all look the same though)
qp_old$lin.fits <- 
  qp_old$summary %>% 
  filter(is.std) %>% 
  group_by(plate,time) %>% # these are all the col variables that distinguish the std curves.
  do(fit = lm(Cq ~ log(conc), data = .) ) %>%  # fit a line to each std curve using lm()
  # now this gets ugly, but basically this is pulling out the y.intercept, slope, and R^2 values
  # by using both the tidy() and glance() functions from broom, and then combining them 
  # together with a left_join
  { left_join(
    tidy(.,fit) %>%   # tidy pulls out the slope and y.intercept, but we have to use spread() to make them separate columns
      select(plate,time,term,estimate) %>% 
      spread(term,estimate) %>% rename(y.intercept=`(Intercept)`, slope=`log(conc)`),
    glance(.,fit) %>% # glance pulls out the r.sqrd values 
      select(plate,time,r.squared, adj.r.squared)
  )
  } %>% 
  ungroup() %>%
  mutate(  qpcr.effic = exp(-1/slope) - 1  ) # qPCR efficiency = e^(-1/slope) - 1

# Use linear fits to calculate concentration (calc.conc) from Cq:
qp_old$summary <- 
  left_join(qp_old$summary, 
            qp_old$lin.fits) %>%
  select(-r.squared,-adj.r.squared,-qpcr.effic) %>% # only keep the added columns for y.intercept and slope
  mutate(calc.conc = case_when( 
    .$is.empty  ~ NA_real_, # empty wells have a calc.conc of NA
    .$is.std    ~ NA_real_, # std curve wells have a calc.conc of NA
    is.na(.$Cq) ~ 0,        # all other wells with Cq=NA have a calc.conc of 0
    TRUE        ~ exp( (.$Cq-.$y.intercept)/.$slope) )) # this is the formula to convert from Cq to calc.conc

# outlier wells (see extra plots for all triplicate measures ... F06 and F07 are
# outliers, likely because I pipetted the DNA that should have gone in F06 into
# F07 by mistake, giving one well 2x the DNA and one well ~0 DNA)
tmp.outlierwells <- c("F06","F07")

# Now do the WT %rad52 loading over time

# This helper function extracts only certin cell cycles that are valid for % Rad52 analysis
# (same as % Rad52 bar plots below, require at least 100 min of usable movie after nuc birth)
df.rad52.valid <- find_pct_rad52_valid_cells(df, min.time = 100)

# This helper fcn takes these valid cell cycles and bins various events, like Rad52 loading
# and lacO dissaperance
# bin.max is the latest usable timepoint across all movies, in minutes:
bin.min <- 0 # minutes
bin.max <- df.rad52.valid %>% mutate(max = (f.last_good_frame-1)*m.t_int + e.t_start) %>% {max(.$max)/60} # minutes
bin.width <- 8 # minutes
rad52.cumulative.counts <- find_cdf_like_rad52_cell_counts(df.rad52.valid, bin.width, bin.min, bin.max)

# combined qPCR AND Rad52% together with two y-axes:
tmp.qp.data <- 
  qp_old$summary %>%
  filter(!is.empty, !is.noDNA, !is.std, primers %in% c("-168bps","Ncb2"), !(well %in% tmp.outlierwells)) %>%
  # replace "-168bps" with "m168" to make variable names easier to handle
  mutate(primers = as.factor(  if_else(primers == "-168bps","m168",as.character(primers))  )) %>%
  select(strain,DNA,time,conc,primers,digest,rep,calc.conc) %>%
  # find % resected for rep 1, 2, and 3 ... then get error bars across reps (some are NA bc of tmp.outlierwells)
  unite(pcr.cond,primers,digest,sep=".") %>%
  spread(pcr.cond,calc.conc) %>%
  mutate(pct.resected    = (m168.ApoI/m168.HincII) - (Ncb2.ApoI/m168.HincII)) %>%
  group_by(time) %>%
  summarize(pct.resected.ave = mean(pct.resected,na.rm=T),
            pct.resected.n   = sum(!is.na(pct.resected)),
            pct.resected.se = sd(pct.resected,na.rm=T) / sqrt(pct.resected.n)) %>%
  mutate(time = if_else(time==110,120,time))

tmp.rad52.cum.data <- 
  rad52.cumulative.counts %>% 
  filter(e.genotype == "wt", event.type == "cmpct_cyc_with_Rad52")

tmp.color.y.right <- "brown4"
tmp.g <- ggplot() +
  geom_bar(data=tmp.qp.data, 
           aes(x=time, y=pct.resected.ave),
           width=46, stat="identity", color="black", fill="gray65", size=0.6) +
  geom_errorbar(data=tmp.qp.data,
                aes(x=time, ymin=pct.resected.ave-pct.resected.se, ymax=pct.resected.ave+pct.resected.se), 
                size = 1, width = 20) + 
  geom_step(data=tmp.rad52.cum.data, 
            aes(x=bin.right, y=g.cm.value),
            size=1.2, color=tmp.color.y.right) + 
  scale_x_continuous(name="Time after HO induction (min)",limits=c(-23,240),breaks=c(0,60,120,180,240,300)) + 
  scale_y_continuous(name="ApoI protection\n(168 nt from HOcs)",limits=c(0,.17),labels=scales::percent,
                     sec.axis = sec_axis(~., name="Cumulative cell cycles\nwith Rad52 foci",labels=scales::percent)) + 
  gg$my.theme +
  theme(axis.title.y.right=element_text(color=tmp.color.y.right),
        axis.text.y.right=element_text(color=tmp.color.y.right))
tmp.g

# use grid graphics to change just the color of the right y-axis!
tmp.grob <- ggplotGrob(tmp.g)
tmp.grob$grobs[[grep("axis-r",tmp.grob$layout$name)]]$children[[1]]$gp$col <- tmp.color.y.right
tmp.grob$grobs[[grep("axis-r",tmp.grob$layout$name)]]$children[[2]][[1]][[1]]$gp$col <- tmp.color.y.right
grid.newpage()
grid.draw(tmp.grob)

ggsave(filename = str_c(fig_save_dir,"Figure1s3C_WT_qPCR_and_PctRad52_overlaid.pdf"),
       plot = grid.draw(tmp.grob),
       device = cairo_pdf(),
       scale = 1, width = 16, height = 11.5, units = "cm")





# qPCR showing resection rate on LacO and non-LacO side of break (Fig1s3H) --------
tmp.to.plot <- 
  qp$summ.calc.conc %>%
  # remove excluded values, and put ApoI and HincII calculated concentrations into two different cols:
  filter(is.expWell, !is.excluded) %>%
  select(date, plate, strain, time, digest, primers, rep, calc.conc) %>%
  spread(digest, calc.conc) %>% 
  # calculate means for each (ApoI, HincII), then get pct.undigested=ApoI/HincII, using error propagation!
  group_by(date, plate, strain, time, primers) %>%
  summarize(ave.ApoI = mean(ApoI,na.rm=T), n.ApoI = sum(!is.na(ApoI)), se.ApoI = sd(ApoI,na.rm=T) / sqrt(n.ApoI),
            ave.HincII = mean(HincII,na.rm=T), n.HincII = sum(!is.na(HincII)), se.HincII = sd(HincII,na.rm=T) / sqrt(n.HincII),
            ave.pct.undigested = ave.ApoI/ave.HincII,
            se.pct.undigested = ave.pct.undigested * sqrt( (se.ApoI/ave.ApoI)^2 + (se.HincII/ave.HincII)^2 )) %>%
  # find the proper Ncb2 controls to use for each measurment
  group_by(date, plate, strain, time) %>%
  mutate(ave.ncb2.digest.eff = first( if_else(primers == "Ncb2_v1",ave.pct.undigested,NA_real_) ),
         se.ncb2.digest.eff = first( if_else(primers == "Ncb2_v1",se.pct.undigested,NA_real_) )) %>%
  filter(primers != "Ncb2_v1") %>% # get rid of all Ncb2 rows, we don't need them anymore!
  # if there is not an Ncb2 control for a specific strain/primer combo, then use an average from the plate instead (error prop. is wonky here ...)
  ungroup() %>% group_by(date,plate) %>%
  mutate(ave.ncb2.digest.eff = if_else(is.finite(ave.ncb2.digest.eff), ave.ncb2.digest.eff, mean(ave.ncb2.digest.eff, na.rm=TRUE)),
         se.ncb2.digest.eff = if_else(is.finite(se.ncb2.digest.eff), se.ncb2.digest.eff, mean(se.ncb2.digest.eff, na.rm=TRUE))) %>%
  # using Ncb2 digest efficiency, "background subtract" and use error propagation
  ungroup() %>%
  mutate(ave.pct.resection = ave.pct.undigested - ave.ncb2.digest.eff,
         se.pct.resection = sqrt(se.pct.undigested^2 + se.ncb2.digest.eff^2))

# with 95% CIs for error bars:
tmp.to.plot %>%
  filter(date == "2017_0727", strain == 2123) %>%
  arrange(primers,time) %>%
  ggplot(aes(x=time, color=primers,
             y=ave.pct.resection, ymin=ave.pct.resection-se.pct.resection*qnorm(0.975), ymax=ave.pct.resection+se.pct.resection*qnorm(0.975))) +
  geom_line(size = 1.25) + geom_point() + 
  geom_errorbar(size = 0.75, width = 10, color="black") + 
  scale_y_continuous(name="Percent Resected\n(undigested by ApoI)", labels=scales::percent) + 
  scale_x_continuous(name="Time after HO induction", breaks = seq(-360,360,120), labels = function(x) abs(x)) +
  scale_color_gdocs(name="Distance from HOcs", labels=c("-168 bps","-14,253 bps","+13,150 bps")) + 
  gg$my.theme

ggsave(filename = str_c(fig_save_dir,"Figure1s3H_qPCR_equal_rates_LacOside_and_other_side.pdf"),
       device = cairo_pdf(),
       scale = 1, width = 20, height = 12.5, units = "cm")











# qPCR at 90min: wt, exo1, ctp1 (Fig 1s3I) --------------------------------
tmp.to.plot <- 
  qp$summ.calc.conc %>%
  # remove excluded values, and put ApoI and HincII calculated concentrations into two different cols:
  filter(is.expWell, !is.excluded) %>%
  select(date, plate, strain, time, digest, primers, rep, calc.conc) %>%
  spread(digest, calc.conc) %>% 
  # calculate means for each (ApoI, HincII), then get pct.undigested=ApoI/HincII, using error propagation!
  group_by(date, plate, strain, time, primers) %>%
  summarize(ave.ApoI = mean(ApoI,na.rm=T), n.ApoI = sum(!is.na(ApoI)), se.ApoI = sd(ApoI,na.rm=T) / sqrt(n.ApoI),
            ave.HincII = mean(HincII,na.rm=T), n.HincII = sum(!is.na(HincII)), se.HincII = sd(HincII,na.rm=T) / sqrt(n.HincII),
            ave.pct.undigested = ave.ApoI/ave.HincII,
            se.pct.undigested = ave.pct.undigested * sqrt( (se.ApoI/ave.ApoI)^2 + (se.HincII/ave.HincII)^2 )) %>%
  # find the proper Ncb2 controls to use for each measurment
  group_by(date, plate, strain, time) %>%
  mutate(ave.ncb2.digest.eff = first( if_else(primers == "Ncb2_v1",ave.pct.undigested,NA_real_) ),
         se.ncb2.digest.eff = first( if_else(primers == "Ncb2_v1",se.pct.undigested,NA_real_) )) %>%
  filter(primers != "Ncb2_v1") %>% # get rid of all Ncb2 rows, we don't need them anymore!
  # if there is not an Ncb2 control for a specific strain/primer combo, then use an average from the plate instead (error prop. is wonky here ...)
  ungroup() %>% group_by(date,plate) %>%
  mutate(ave.ncb2.digest.eff = if_else(is.finite(ave.ncb2.digest.eff), ave.ncb2.digest.eff, mean(ave.ncb2.digest.eff, na.rm=TRUE)),
         se.ncb2.digest.eff = if_else(is.finite(se.ncb2.digest.eff), se.ncb2.digest.eff, mean(se.ncb2.digest.eff, na.rm=TRUE))) %>%
  # using Ncb2 digest efficiency, "background subtract" and use error propagation
  ungroup() %>%
  mutate(ave.pct.resection = ave.pct.undigested - ave.ncb2.digest.eff,
         se.pct.resection = sqrt(se.pct.undigested^2 + se.ncb2.digest.eff^2))

tmp.dist.factor.levels <- c("m300","m3023") # "m1556" ?
tmp.dist.factor.labels <- c("-300 bps", "-3,023 bps")

tmp.geno.factor.levels <- c(2123,1914,2477)
tmp.geno.factor.labels <- c("paste('WT')","paste(italic('exo1'),Delta)","paste(italic('ctp1'),Delta)")


# tmp.to.plot is calculated %resection from above "New 1E"
tmp.to.plot %>%
  filter(strain %in% tmp.geno.factor.levels,
         primers %in% tmp.dist.factor.levels,
         time %in% c(0,90)) %>%
  arrange(strain,primers,time) %>%
  mutate(strain = factor(strain, tmp.geno.factor.levels, tmp.geno.factor.labels),
         primers = factor(primers, tmp.dist.factor.levels, tmp.dist.factor.labels)) %>%
  
  # normalize to set 0min %resection to zero:
  # this is dangerous bc it remove the exp data/plate info ...
  select(strain,time,primers,  ave.pct.resection,se.pct.resection) %>% 
  gather(key = "ave_or_se", value = "value", ave.pct.resection, se.pct.resection) %>%
  unite(new_key,ave_or_se,time,sep="_") %>%
  spread(new_key,value) %>%
  
  # now calculate the new normalized pct resection at 90 min:
  mutate(norm.90min.pct.resection = ave.pct.resection_90 - ave.pct.resection_0,
         norm.90min.se = sqrt( se.pct.resection_0^2 + se.pct.resection_90^2 )) %>%
  
  # plot
  ggplot(aes(x=primers, 
             y=norm.90min.pct.resection, ymin=norm.90min.pct.resection-norm.90min.se*qnorm(0.975), ymax=norm.90min.pct.resection+norm.90min.se*qnorm(0.975))) +
  geom_bar(stat="identity",position="dodge",color="black", fill="grey") + 
  geom_errorbar(size = 1, width = 0.5) + 
  geom_hline(yintercept = 0) + 
  facet_grid(.~strain, scales="free", space="free", labeller = label_parsed) + 
  scale_y_continuous(name="Percent Resected\n(normalized to T=0min)", labels=scales::percent) + 
  scale_x_discrete(name="Distance from HOcs") + 
  gg$my.theme + 
  gg$rotate.x.lab + 
  theme(strip.text = element_text(size = 20, face = "bold"))

ggsave(filename = str_c(fig_save_dir,"Figure1s3I_qPCR_wt_exo1_ctp1_90normalized.pdf"),
       device = cairo_pdf(),
       scale = 1, width = 12, height = 13, units = "cm")



# Pct Cells with Rad52 foci at 180min (Fig1s3J) ------------------------------


# valid cells based on all critera above using this helper function
df.rad52.valid <- find_pct_rad52_valid_cells(df, min.time = 100)

# use this helper function to calculate the means and error bars for % Rad52 on
# both a per-experiment (biological replicate) basis (columns starting with "e.") 
# as well as on a per-genotype basis ("g." columns). Includes weighted, pooled, 
# and simple averaged means and variance.
min.cells.cutoff <- 35
pct.rad52.by.run <- calc_rad52_pcts(df.rad52.valid, min.cells.cutoff = min.cells.cutoff)

# genotypes we are interested in for this figure:
tmp.genotypes <- c("wt","exo1")

# All data pooled across biological replicates, 95% CI error bars:
pct.rad52.by.run %>%
  filter(e.genotype %in% tmp.genotypes) %>%
  ggplot(aes(x=e.genotype, y=g.pooled_mean, ymin=g.pooled_mean-g.pooled_ci95, ymax=g.pooled_mean+g.pooled_ci95)) +
  geom_bar(stat="unique",position="identity") + 
  geom_errorbar(size = 1, width = 0.35) + 
  scale_x_discrete(name="",labels=gg$genotype.labels,drop=TRUE) +
  scale_y_continuous(name="Cell cycles w/ Rad52 foci",limits=c(0,.2),labels=scales::percent) + 
  gg$my.theme + gg$rotate.x.lab

ggsave(filename = str_c(fig_save_dir,"Figure1s3J_Rad52_foci_formation_WT_exo1.pdf"),
       device = cairo_pdf(),
       scale = 1, width = 7, height = 12, units = "cm")




# qPCR comparison of WT and rev (Fig 2D) ----------------------------------

tmp.to.plot <- 
  qp$summ.calc.conc %>%
  # remove excluded values, and put ApoI and HincII calculated concentrations into two different cols:
  filter(is.expWell, !is.excluded) %>%
  select(date, plate, strain, time, digest, primers, rep, calc.conc) %>%
  spread(digest, calc.conc) %>% 
  # calculate means for each (ApoI, HincII), then get pct.undigested=ApoI/HincII, using error propagation!
  group_by(date, plate, strain, time, primers) %>%
  summarize(ave.ApoI = mean(ApoI,na.rm=T), n.ApoI = sum(!is.na(ApoI)), se.ApoI = sd(ApoI,na.rm=T) / sqrt(n.ApoI),
            ave.HincII = mean(HincII,na.rm=T), n.HincII = sum(!is.na(HincII)), se.HincII = sd(HincII,na.rm=T) / sqrt(n.HincII),
            ave.pct.undigested = ave.ApoI/ave.HincII,
            se.pct.undigested = ave.pct.undigested * sqrt( (se.ApoI/ave.ApoI)^2 + (se.HincII/ave.HincII)^2 )) %>%
  # find the proper Ncb2 controls to use for each measurment
  group_by(date, plate, strain, time) %>%
  mutate(ave.ncb2.digest.eff = first( if_else(primers == "Ncb2_v1",ave.pct.undigested,NA_real_) ),
         se.ncb2.digest.eff = first( if_else(primers == "Ncb2_v1",se.pct.undigested,NA_real_) )) %>%
  filter(primers != "Ncb2_v1") %>% # get rid of all Ncb2 rows, we don't need them anymore!
  # if there is not an Ncb2 control for a specific strain/primer combo, then use an average from the plate instead (error prop. is wonky here ...)
  ungroup() %>% group_by(date,plate) %>%
  mutate(ave.ncb2.digest.eff = if_else(is.finite(ave.ncb2.digest.eff), ave.ncb2.digest.eff, mean(ave.ncb2.digest.eff, na.rm=TRUE)),
         se.ncb2.digest.eff = if_else(is.finite(se.ncb2.digest.eff), se.ncb2.digest.eff, mean(se.ncb2.digest.eff, na.rm=TRUE))) %>%
  # using Ncb2 digest efficiency, "background subtract" and use error propagation
  ungroup() %>%
  mutate(ave.pct.resection = ave.pct.undigested - ave.ncb2.digest.eff,
         se.pct.resection = sqrt(se.pct.undigested^2 + se.ncb2.digest.eff^2))


tmp.dist.factor.levels <- c("m300","m3023","p13150") # "m1556" ?
tmp.dist.factor.labels <- c("-300 bps", "-3,023 bps", "+13,150 bps")

tmp.geno.factor.levels <- c(2123,2149)
tmp.geno.factor.labels <- c("wt","rev7")
tmp.90min.from.Fig4C <- 
  tmp.to.plot %>%
  filter(!(date == "2017_0807" & strain == 2149),
         strain %in% tmp.geno.factor.levels,
         primers %in% c("m300","m3023"),
         time %in% c(0,90)) %>%
  arrange(strain,primers,time) %>%
  mutate(strain = factor(strain, tmp.geno.factor.levels, tmp.geno.factor.labels),
         primers = factor(primers, tmp.dist.factor.levels, tmp.dist.factor.labels)) %>% 
  
  # normalize to set 0min %resection to zero:
  # this is dangerous bc it remove the exp data/plate info ...
  select(strain,time,primers,  ave.pct.resection,se.pct.resection) %>% 
  gather(key = "ave_or_se", value = "value", ave.pct.resection, se.pct.resection) %>%
  unite(new_key,ave_or_se,time,sep="_") %>%
  spread(new_key,value) %>%
  
  # now calculate the new normalized pct resection at 90 min:
  mutate(time.facet.label = "090min",
         norm.pct.resection = ave.pct.resection_90 - ave.pct.resection_0,
         norm.se = sqrt( se.pct.resection_0^2 + se.pct.resection_90^2 ))


tmp.180min.from.Fig2C <- 
  tmp.to.plot %>%
  filter(date == "2017_0807", strain %in% c(2123,2149), primers == "p13150") %>%
  arrange(strain,primers,time) %>%
  mutate(strain = factor(strain, tmp.geno.factor.levels, tmp.geno.factor.labels),
         primers = factor(primers, tmp.dist.factor.levels, tmp.dist.factor.labels)) %>% 
  
  # normalize to set 0min %resection to zero:
  # this is dangerous bc it remove the exp data/plate info ...
  select(strain,time,primers,  ave.pct.resection,se.pct.resection) %>% 
  gather(key = "ave_or_se", value = "value", ave.pct.resection, se.pct.resection) %>%
  unite(new_key,ave_or_se,time,sep="_") %>%
  spread(new_key,value) %>%
  
  # now calculate the new normalized pct resection at 90 min:
  mutate(time.facet.label = "180min",
         norm.pct.resection = ave.pct.resection_180 - ave.pct.resection_0,
         norm.se = sqrt( se.pct.resection_0^2 + se.pct.resection_180^2 ))

tmp.comb.to.plot <- 
  bind_rows(tmp.90min.from.Fig4C  %>% select(-ave.pct.resection_0,-ave.pct.resection_90,-se.pct.resection_0,-se.pct.resection_90),
            tmp.180min.from.Fig2C %>% select(-ave.pct.resection_0,-ave.pct.resection_180,se.pct.resection_0,-se.pct.resection_180))


tmp.comb.to.plot %>%
  mutate(facet.var = interaction(primers,time.facet.label)) %>%
  ggplot(aes(x=strain, fill = strain,
             y=norm.pct.resection, ymin=norm.pct.resection-norm.se*qnorm(0.975), ymax=norm.pct.resection+norm.se*qnorm(0.975))) +
  geom_bar(stat="identity",position="dodge",color="black") + 
  geom_errorbar(size = 1, width = 0.5) + 
  geom_hline(yintercept = 0) + 
  facet_grid(.~facet.var, scales="free", space="free") + 
  scale_y_continuous(name="Percent Resected\n(normalized to T=0min)", labels=scales::percent) + 
  scale_x_discrete(name="") + 
  scale_fill_gdocs(name="", labels=c("WT",expression(paste(italic("rev7"),Delta)))) +
  coord_cartesian(xlim = c(0.25,2.75), expand = FALSE) + 
  gg$my.theme + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank())

ggsave(filename = str_c(fig_save_dir,"Figure2D_qPCR_wt_and_rev7_comparison.pdf"),
       device = cairo_pdf(),
       scale = 1, width = 18, height = 12, units = "cm")
