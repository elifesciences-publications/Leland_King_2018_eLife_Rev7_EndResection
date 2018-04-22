
# Imports and packages ----------------------------------------------------
library(tidyverse)
library(stringr)
library(readxl)
library(broom)
library(ggthemes)


# Function for reading in qPCR data ---------------------------------------

# Read in .xlsx files with all raw data from one or more Bio-Rad qPCR runs
# and make a list called "qp" which has just a few tidy data frames: run.info, 
# summary, cycle, and melt
qpcr_import <- function(...) {
  
  # validate an even number of inputs:
  inp_args <- list(...)
  if ( length(inp_args) == 0 || 
       (length(inp_args) %% 2) != 0 ) {
    stop("Number of inputs must be even. Give inputs in pairs: xlsx.dir1, well.info.file1, xlsx.dir2, well.info.file2, ...")
  }
  # convert to more readable table:
  inp_args <- inp_args %>% unlist() %>% matrix(nrow=2,dimnames=list(c("xlsx.dir","well.info.file"))) %>% t() %>% as_tibble()
  
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
    well.info <- read_excel(str_c(xlsx.dir,well.info.file),sheet = 1)
    
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
    qp$run.info <- read_excel(str_c(xlsx.dir,files.qpcr["quant.sum"]), sheet=2, col_names=FALSE)
    
    
    # build up the qp$summary df, which will have only one row per well, so it will contain
    # all data that is one per well (like Cq, end point, Melt Peak, etc.)
    tmp.cq <- read_excel(str_c(xlsx.dir,files.qpcr["quant.sum"]),sheet=1) # get well, fluor and Cq values from "Quantification Summary"
    tmp.cq <- tmp.cq %>% select(well=Well,fluor=Fluor,Cq) %>% 
      mutate(well = as.factor(letter_numb_pad_zeros(well)),
             fluor = as.factor(fluor))
    
    tmp.endpts <- read_excel(str_c(xlsx.dir,files.qpcr["end.points"]),sheet=1) # get end.RFU from "End Point Results"
    tmp.endpts <- tmp.endpts %>% select(well=Well,end.RFU=`End RFU`) %>% 
      mutate(well = as.factor(letter_numb_pad_zeros(well)))
    
    tmp.melt <- read_excel(str_c(xlsx.dir,files.qpcr["melt.peaks"]),sheet=1) # get melt-related info from "Melt Curve Peak Results"
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
    tmp.cycle <- read_excel(str_c(xlsx.dir,files.qpcr["quant.amp"]),sheet=1)
    qp$cycle <- 
      tmp.cycle %>%
      gather(well,RFU,-Cycle) %>%
      mutate(well = as.factor(letter_numb_pad_zeros(well))) %>%
      select(well, RFU, cycle = Cycle) %>%
      left_join(qp$summary,by="well") %>%
      select(everything(),-cycle,-RFU,cycle,RFU)
    
    # Make qp$melt (has FULL traces for each melt curve ... many rows of data per well!)
    tmp.meltamp <- read_excel(str_c(xlsx.dir,files.qpcr["melt.amp"]),sheet=1) # ger raw melt curve data (RFU value at each temp for each well)
    tmp.meltamp <- tmp.meltamp %>%
      gather(well,RFU,-Temperature) %>%
      mutate(well = as.factor(letter_numb_pad_zeros(well))) %>%
      select(well, temp = Temperature, RFU)
    
    tmp.meltder <- read_excel(str_c(xlsx.dir,files.qpcr["melt.der"]),sheet=1) # get derrivative of melt curve
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
    
    # add a plate column (=r) to left of all data frames:
    this.qp$run.info <- this.qp$run.info %>% mutate(plate = as.integer(r)) %>% select(plate,field=1,value=2)
    this.qp$summary  <- this.qp$summary  %>% mutate(plate = as.integer(r)) %>% select(plate,everything())
    this.qp$cycle    <- this.qp$cycle    %>% mutate(plate = as.integer(r)) %>% select(plate,everything())
    this.qp$melt     <- this.qp$melt     %>% mutate(plate = as.integer(r)) %>% select(plate,everything())
    
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



# Import the raw qPCR data ------------------------------------------------

tmp.rdata.file.name <- "2017_0706_all_raw_data.RData"
if ( file.exists(tmp.rdata.file.name) ) {
  load(tmp.rdata.file.name)
} else {
  # use the helper function qpcr_import() from above to import the raw excel data:
  qp <- qpcr_import("RawData/","2017_0706_well_info_MCK.xlsx")
  # save this as an .RData file to speed things up next time around:
  save(qp, file = tmp.rdata.file.name )
}



# Visualize some sanity checks first: -------------------------------------

# std curves with ALL points included
tmp.base <-  10
qp$summary %>%
  filter(is.std) %>%
  mutate(uniq.std.curves = interaction(primers,digest,sep=" ")) %>%
  ggplot(aes(x=log(conc,tmp.base), y=Cq, group=uniq.std.curves, color=uniq.std.curves)) + 
  geom_point(size=2, stroke=2, shape=21) + 
  facet_wrap(~uniq.std.curves) + 
  stat_smooth(geom = "line",method = "lm") + 
  scale_color_few(name="Std Curve") + 
  scale_y_continuous(name="Cq value") +
  scale_x_continuous(name=str_c("Concentration [ng/uL], log",tmp.base)) + 
  theme_bw()
ggsave(filename = file.path("Rplots","validate_std_curves.pdf"),
       device = cairo_pdf(),
       scale = 1, width = 14, height = 10, units = "cm")

qp$cycle %>%
  filter(is.std) %>%
  mutate(uniq.std.curves = interaction(primers,digest,sep=" ")) %>%
  ggplot(aes(x=cycle, y=RFU, group=well, color=uniq.std.curves)) + 
  geom_line(size=0.4) +
  facet_wrap(~uniq.std.curves) + 
  scale_color_few(name="Std Curve") + 
  scale_y_continuous(name="RFU") +
  scale_x_continuous(name="Cycle") + 
  theme_bw()
ggsave(filename = file.path("Rplots","validate_std_curve_rfu_plots.pdf"),
       device = cairo_pdf(),
       scale = 1, width = 14, height = 10, units = "cm")



# all amplifications curves (all wells are non-std curve wells for this exp):
qp$cycle %>%
  ggplot(aes(x=cycle, y=RFU, group=well)) + 
  geom_line(size=0.2) +
  facet_wrap(~primers) + 
  scale_y_continuous(name="RFU") +
  scale_x_continuous(name="Cycle") + 
  theme_bw()
ggsave(filename = file.path("Rplots","validate_check_all_amp_curves.pdf"),
       device = cairo_pdf(),
       scale = 1, width = 14, height = 14, units = "cm")

# all melt curves for all +DNA wells:
qp$melt %>%
  filter(!is.noDNA) %>%
  ggplot(aes(x=temp, y=dRFU, group=well)) + 
  geom_line(size=0.4, color="grey75") +
  facet_wrap(~primers) + 
  scale_y_continuous(name="RFU") +
  scale_x_continuous(name="Temperature (C)") + 
  theme_bw()
ggsave(filename = file.path("Rplots","validate_check_primer_melt_peaks.pdf"),
       device = cairo_pdf(),
       scale = 1, width = 14, height = 10, units = "cm")



# Std curves and fits: ----------------------------------------------------

qp$lin.fits <- 
  qp$summary %>% 
  filter(is.std) %>% 
  group_by(primers,digest) %>% # these are all the col variables that distinguish the std curves.
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

# view and save the linear fits:
qp$lin.fits
write_csv(qp$lin.fits, file.path("Rplots","std_curve_fit_stats.csv"))

# Use linear fits to calculate concentration (calc.conc) from Cq:
qp$summ.calc.conc <- 
  qp$lin.fits %>%
  filter(digest == "HincII") %>% # need to get rid of the "none" std curve and digest column before joining!!!
  select(-digest) %>%
  left_join(qp$summary,.) %>% 
  mutate(calc.conc = case_when(.$is.empty  ~ NA_real_, # empty wells have a calc.conc of NA
                               is.na(.$Cq) ~ 0,        # all other wells with Cq=NA have a calc.conc of 0
                               TRUE        ~ exp( (.$Cq-.$y.intercept)/.$slope) )) # this is the formula to convert from Cq to calc.conc w/ natural log



# Bar plots!! -------------------------------------------------------------

# All Cq values for every relavent well
tmp.exp.conc <- 2.5
qp$summary %>%
  filter(conc == tmp.exp.conc) %>% # get rid of everything except the ApoI protection exp
  arrange(strain,primers,time,digest) %>% 
  mutate(bar.groups = str_c(strain," ",primers," ",time,"min"," ",digest)) %>% 
  ggplot(aes(x=bar.groups, y=Cq, fill=factor(rep))) +
  geom_bar(stat="identity",position="dodge",color="black") + 
  scale_fill_manual(values = c("grey","grey","grey"),guide=FALSE) + 
  scale_y_continuous(name="Cq value") + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = rel(1), angle = 45, hjust = 1, color = "black"))

# Exclude these wells:
#qp$summ.calc.conc %>% 
#  filter(primers == "m3023", digest == "HincII", time == 0, conc == 2.5) %>%
#  select(plate:is.empty,Cq,calc.conc)
tmp.excluded.wells <- c("C05")

# All Cq values for every relavent well (mark excluded ones though):
qp$summary %>%
  filter(conc == tmp.exp.conc) %>% # get rid of everything except the ApoI protection exp
  mutate(is.excluded = if_else(well %in% tmp.excluded.wells, TRUE, FALSE)) %>%
  arrange(strain,primers,time,digest,rep) %>% 
  mutate(bar.groups = str_c(strain," ",primers," ",time,"min"," ",digest," ",rep)) %>% 
  ggplot(aes(x=bar.groups, y=Cq, fill=is.excluded)) +
  geom_bar(stat="identity",color="black") + 
  scale_fill_manual(name="excluded wells",values = c("grey",NA)) + 
  scale_y_continuous(name="Cq value") + 
  scale_x_discrete(name="") + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = rel(1), angle = 90, hjust = 1, color = "black"))

ggsave(filename = file.path("Rplots","bar_plot_all_raw_Cq_values.pdf"),
       device = cairo_pdf(),
       scale = 1, width = 22, height = 14, units = "cm")


# Cq means + error bars:
qp$summary %>%
  filter(conc == tmp.exp.conc) %>% # get rid of everything except the ApoI protection exp
  mutate(is.excluded = if_else(well %in% tmp.excluded.wells, TRUE, FALSE)) %>%
  arrange(strain,primers,time,digest,rep) %>% 
  mutate(bar.groups = str_c(strain," ",primers," ",time,"min"," ",digest)) %>% 
  group_by(bar.groups) %>%
  summarize(mean.Cq = mean(Cq,na.rm=T), sd = sd(Cq,na.rm=T), n = sum(!is.na(Cq)), se = sd/sqrt(n)) %>%
  ggplot(aes(x=bar.groups, y=mean.Cq, ymin=mean.Cq-se, ymax=mean.Cq+se)) +
  geom_bar(stat="identity",color="black", fill="grey") + 
  geom_errorbar(size = 0.2, width = 0.5) + 
  scale_y_continuous(name="mean Cq value") + 
  scale_x_discrete(name="") + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = rel(1), angle = 90, hjust = 1, color = "black"))



# Calculate % resected using ApoI protection
# NOTE: this does NOT have any std curves, so this assumes 100%
# qPCR efficiency and perfectly linearly Cq values for everything!
tmp.exp.conc <- 2.5
qp$summary %>%
  filter(conc == tmp.exp.conc) %>% # get rid of everything except the ApoI protection exp
  mutate(is.excluded = if_else(well %in% tmp.excluded.wells, TRUE, FALSE)) %>%
  select(strain,time,primers,rep,digest,Cq) %>%
  spread(digest,Cq) %>%
  mutate(pct.undigested = 2^-(ApoI - HincII)) %>%
  group_by(strain,time) %>%
  mutate(ncb2.digest.eff = mean( if_else(primers == "Ncb2_v1",pct.undigested,NA_real_) ,na.rm=TRUE)) %>%
  filter(primers != "Ncb2_v1") %>%
  mutate(pct.resection = pct.undigested - ncb2.digest.eff) %>%
  group_by(strain,time,primers) %>%
  summarize(mean.pct.resection = mean(pct.resection,na.rm=T), 
            sd.pct.resection= sd(pct.resection,na.rm=T),
            n = sum(!is.na(pct.resection)),
            serr.pct.resection = sd.pct.resection/sqrt(n)) %>%
  arrange(strain,time,primers) %>% 
  mutate(bar.groups = str_c(strain," ",time,"min"," ",primers)) %>%
  
  ggplot(aes(x=bar.groups, 
             y=mean.pct.resection,
             ymin=mean.pct.resection-serr.pct.resection,
             ymax=mean.pct.resection+serr.pct.resection)) +
  geom_bar(stat="identity",position="dodge",color="black", fill="grey") + 
  geom_errorbar(size = 1, width = 0.5) + 
  scale_y_continuous(name="Percent Resected", labels=scales::percent) + 
  scale_x_discrete(name="") + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = rel(1), angle = 45, hjust = 1, color = "black"))

ggsave(filename = file.path("Rplots","bar_plot_pct_resection_NoStdCurve.pdf"),
       device = cairo_pdf(),
       scale = 1, width = 14, height = 12, units = "cm")




# TESTING
tmp.exp.conc <- 2.5
qp$summ.calc.conc %>%
  filter(conc == tmp.exp.conc) %>%
  mutate(is.excluded = if_else(well %in% tmp.excluded.wells, TRUE, FALSE),
         plot.conc = if_else(is.finite(calc.conc),calc.conc,conc)) %>%
  filter(!is.excluded) %>%
  select(strain, time, primers, digest, rep, plot.conc) %>%
  spread(digest, plot.conc) %>%
  mutate(pct.undigested = ApoI/HincII) %>%
  group_by(strain,time) %>%
  mutate(ncb2.digest.eff = mean( if_else(primers == "Ncb2_v1",pct.undigested,NA_real_) ,na.rm=TRUE),
         pct.resection = pct.undigested - ncb2.digest.eff) %>%
  filter(primers != "Ncb2_v1") %>%
  group_by(strain, time, primers) %>%
  summarize(mean.pct.resection = mean(pct.resection,na.rm=T), 
            sd.pct.resection = sd(pct.resection,na.rm=T),
            n = sum(!is.na(pct.resection)),
            serr.pct.resection = sd.pct.resection / sqrt(n) ) %>%
  arrange(strain,primers,time) %>% 
  mutate(bar.groups = str_c(strain," ",primers," ",time,"min")) 
