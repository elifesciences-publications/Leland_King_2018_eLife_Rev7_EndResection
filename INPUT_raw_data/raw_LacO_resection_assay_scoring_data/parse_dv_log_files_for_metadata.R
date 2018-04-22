require(tidyverse)
require(readxl)
require(stringr)
require(stringi)
require(lubridate)

# Generate a tibble with metadata for each field

# Some of these variables will be constant throughout one experiment
# but some of them (like xy position) are unique per filed
# Ill go through and confirm that each experiment is identical b/t
# fields for the "m. ..." variables that should be the same
# for all fields of one exp



# Read in "fields" tab of the excel sheet ---------------------------------

xl <- "../raw_scoring_data/2016_blinded_analysis.xlsx"

# fields:
flds <- read_excel(xl, sheet = "f_randorder_aug16")
# cannot have ANY missing data here!
stopifnot( all(!is.na(flds)) )
# must all be unique combinations!
stopifnot( nrow(flds) == nrow(distinct(flds)) )


# Add new columns to this tibble ------------------------------------------

# these will be populated by the for loop below
flds <- flds %>% mutate(m.t_num     = as.numeric(NA),    # check that t_num*3*z_num = total number of images acquired!
                        m.t_int     = as.duration(NA),    # check that t_int*(t_num-1) = total delta time
                        m.t_tot     = as.duration(NA),
                        
                        m.z_num     = as.numeric(NA),
                        m.z_size    = as.numeric(NA),
                        
                        m.binning   = as.numeric(NA),
                        m.camera    = as.character(NA),
                        m.gain      = as.numeric(NA),
                        
                        m.autofocus = as.character(NA),  # can use tidyr separate() to get this into separate columns 
                        m.fastacq   = as.character(NA),  # ditto
                        
                        m.pol_trans = as.numeric(NA),
                        m.pol_exp   = as.numeric(NA),
                        m.gfp_trans = as.numeric(NA),
                        m.gfp_exp   = as.numeric(NA),
                        m.mch_trans = as.numeric(NA),
                        m.mch_exp   = as.numeric(NA),
                        
                        m.stage_x   = as.numeric(NA),
                        m.stage_y   = as.numeric(NA),
                        m.stage_z   = as.numeric(NA))


# Function to parse a single .log file ------------------------------------
folder <- "../../PRJs_only"

for (r in 1:nrow(flds)) {
  this.prj.file <- str_c(folder, flds$e.folder[r], flds$f.file_prj[r],sep="/")
  this.log.file <- str_replace(this.prj.file, "_[A-Z]{3}.dv", ".dv.log")
  
  
  # read in the FULL log file
  file.success <- tryCatch({
    file.connect <- file(this.log.file,"r")
    log <- readLines(file.connect)
    close(file.connect)
    file.success <- TRUE
  },
    error = function(e) FALSE
  )
  if (!file.success) {
    print(str_c("error at ",r))
    next
    }
  

  # split the log into two more manageable parts (for speed/simplicty):
  idx.exp_record <- log %>% str_detect("Experiment Record:") %>% which()
  stopifnot(  length(idx.exp_record)==1  )
  log.macro      <- log[1:idx.exp_record-1]
  log.exp_record <- log[-c(1:idx.exp_record-1)]
  

  
  
  # time dimension
  ptrn <- "^\\s.*ZWT Dimensions \\(expected\\):\\t(\\d+) x (\\d+) x (\\d+)$"
  dims.zwt <- log.macro %>% str_subset(ptrn) %>% str_match(ptrn) %>% .[1,-1] %>% as.numeric()
  flds$m.t_num[r] <- dims.zwt[3]
  
  ptrn <- "^#KEY TOTAL_TIME (\\d+),(\\d+),(\\d+),(\\d+)$"  # given as: h,m,s,ms
  hmsms <- log.macro %>% str_subset(ptrn) %>% str_match(ptrn) %>% .[1,-1] %>% as.numeric()
  flds$m.t_tot[r] <- duration(hours=hmsms[1], minutes=hmsms[2], seconds=hmsms[3] + hmsms[4]/1000)
  
  ptrn <- "^#KEY TIME_LAPSE (\\d+),(\\d+),(\\d+),(\\d+)$"  # given as: h,m,s,ms
  hmsms <- log.macro %>% str_subset(ptrn) %>% str_match(ptrn) %>% .[1,-1] %>% as.numeric()
  flds$m.t_int[r] <- duration(hours=hmsms[1], minutes=hmsms[2], seconds=hmsms[3] + hmsms[4]/1000)


  # TODO, grab the last image numb ... for validation ... drop this later though!
  
  
  # z dimension
  flds$m.z_num[r] <- dims.zwt[1]
  
  ptrn <- "^#KEY SECT_SPACING ([\\d\\.]+)$"
  flds$m.z_size[r] <- log.macro %>% str_subset(ptrn) %>% str_match(ptrn) %>% .[1,-1] %>% as.numeric()
  
  
  # camera
  ptrn <- "^\\s{2}Type:\\s{8}\\t(.+)$"
  flds$m.camera[r] <- log.macro %>% str_subset(ptrn) %>% str_match(ptrn) %>% .[1,-1]
  
  ptrn <- "^\\s{2}Gain:\\s{8}\\t([\\d\\.]+)X$"
  flds$m.gain[r] <- log.macro %>% str_subset(ptrn) %>% str_match(ptrn) %>% .[1,-1] %>% as.numeric()
  
  ptrn <- "^\\s{2}Binning:\\s{18}\\t(\\d+)x\\1$"
  flds$m.binning[r] <- log.macro %>% str_subset(ptrn) %>% str_match(ptrn) %>% .[1,-1] %>% as.numeric()
  
  
  # focusing
  ptrn <- "^\\s{2}FIRSTPOINT AF,([\\d\\.,]+)$"
  flds$m.autofocus[r] <- log.macro %>% str_subset(ptrn) %>% str_match(ptrn) %>% .[1,-1] # note: can use tidyr separate()/unite() on this column!

  ptrn <- "^#KEY DO_FAST_ACQUISITION (T|F)$"
  is.fast_acq <- log.macro %>% str_subset(ptrn) %>% str_match(ptrn) %>% .[1,-1] %>% as.logical()
  if (is.fast_acq) {
    ptrn <- c("^#KEY FAST_TYPE (\\w+)$","^#KEY FAST_SHUTTER_EACH (\\w+)$","^#KEY FAST_READOUT_MODE (\\w+)$","^#KEY LAMPS_OFF_AT_END (\\w+)$")
    # can use separate() if needed on this column. here is the format:
    # fast_acq_type (ex ZWT), fast_shutter_on_each (ex: EXPOSRUE), readout_mode (ex: NORMAL), lamps_off_at_end (ex: )
    flds$m.fastacq[r] <- log.macro %>% str_subset(ptrn) %>% str_match(ptrn) %>% .[,2] %>% str_c(collapse=",")
    
  } else {
    flds$m.fastacq[r] <- as.character(NA)
  }
  
  
  # channels and exposures
  ptrn <- "^#KEY CHANNEL ([\\d\\.]+),(\\w+),\\w+,(\\d+)%"   # some of these end with ",1" but other don't (change in .log format?)
  ch.info_matrix   <- log.macro %>% str_subset(ptrn) %>% str_match(ptrn) %>% .[,c(3,2,4)]
  flds$m.pol_exp[r]   <- dseconds( ch.info_matrix[ch.info_matrix[,1] == "POL",2] )
  flds$m.pol_trans[r] <- as.numeric(ch.info_matrix[ch.info_matrix[,1] == "POL",3]) / 100
  flds$m.gfp_exp[r]   <- dseconds( ch.info_matrix[ch.info_matrix[,1] == "GFP",2] )
  flds$m.gfp_trans[r] <- as.numeric(ch.info_matrix[ch.info_matrix[,1] == "GFP",3]) / 100
  flds$m.mch_exp[r]   <- dseconds( ch.info_matrix[ch.info_matrix[,1] == "mCherry",2] )
  flds$m.mch_trans[r] <- as.numeric(ch.info_matrix[ch.info_matrix[,1] == "mCherry",3]) / 100

  
  # stage position in x,y,z
  ptrn <- "^\\s{11}Stage coordinates:\\s{3}\\(([-+]?[0-9]*\\.?[0-9]+),([-+]?[0-9]*\\.?[0-9]+),([-+]?[0-9]*\\.?[0-9]+)\\)$"
  all.xyz.cords <- log.exp_record %>% str_subset(ptrn) %>% str_match(ptrn) %>% .[,-1] %>% apply(c(1,2),as.numeric)
  # just average each column for now ... could get fancier with this in the future!
  flds$m.stage_x[r] <- mean(all.xyz.cords[,1])
  flds$m.stage_y[r] <- mean(all.xyz.cords[,2])
  flds$m.stage_z[r] <- mean(all.xyz.cords[,3])
  
}



# TODO: validation and fill in missing rows here!!! -----------------------


# TODO: fill in missing rows manually by looking at Bio-Formats metadata!


# TODO: validation here!!!
# drop col with last image number now
flds %>% sapply(FUN=class)

flds %>% select(everything(),-starts_with("f."),-contains("stage")) %>% lapply(FUN = unique)



# Save flds as .RData -----------------------------------------------------

save(flds,file="../../rdata_files/flds_with_metadata.RData")

