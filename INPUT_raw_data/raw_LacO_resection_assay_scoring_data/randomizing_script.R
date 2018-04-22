require(readxl)
require(dplyr)
require(stringr)
require(stringi)



# add in missing metadata to flds ------------------------------------

#Import current randomization work done ~Aug 2016

load("../rdata_files/flds_with_metadata.RData")

# ones that do not have metadata:
flds %>% filter(!complete.cases(.)) %>% select(f.order_aug16,e.folder)

# all fields from these two runs which have some metadata missing:
tmp <- flds %>% filter(e.folder %in% c("2015_0330_run2_1771wt_15fields_60min_day4_7hrAdeNPG","2015_0821_run1_2149rev7"))

# replace all the missing metadata (Excpet stage position) for the two Rev7 fields with missing data:
tmp_idx_cols <- grepl("m.",names(flds)) & !grepl("stage",names(flds))
flds[flds$e.folder == "2015_0821_run1_2149rev7" & !complete.cases(flds), tmp_idx_cols] <- 
  flds[flds$e.folder == "2015_0821_run1_2149rev7" & complete.cases(flds), tmp_idx_cols][1,]

# check again:
tmp <- flds %>% filter(e.folder %in% c("2015_0330_run2_1771wt_15fields_60min_day4_7hrAdeNPG","2015_0821_run1_2149rev7"))

# now manually addin metadata for wt 3/30/15 exp, which has tons of totally corrupted files =(
tmp_idx_rows <- flds$e.folder == "2015_0330_run2_1771wt_15fields_60min_day4_7hrAdeNPG" & !complete.cases(flds)
flds[tmp_idx_rows,"m.t_num"] <- 49
require(lubridate)
flds[tmp_idx_rows,"m.t_int"]     <- duration(minutes=10)
flds[tmp_idx_rows,"m.t_tot"]     <- duration(hours=8)
flds[tmp_idx_rows,"m.z_num"]     <- 25 # this is just a guess, since all I've got is the prj right now ...
flds[tmp_idx_rows,"m.z_size"]    <- 0.26
flds[tmp_idx_rows,"m.binning"]   <- 1
flds[tmp_idx_rows,"m.camera"]    <- "PM1394Cam00 EMCCD / EVOLVE- 512"
flds[tmp_idx_rows,"m.gain"]      <- NA # unknown
flds[tmp_idx_rows,"m.autofocus"] <- NA # unknown
flds[tmp_idx_rows,"m.fastacq"]   <- "ZWT,EXPOSURE,NORMAL,F"
flds[tmp_idx_rows,"m.pol_trans"] <- NA # unknown 
flds[tmp_idx_rows,"m.pol_exp"]   <- NA # unknown 
flds[tmp_idx_rows,"m.gfp_trans"] <- NA # unknown 
flds[tmp_idx_rows,"m.gfp_exp"]   <- NA # unknown 
flds[tmp_idx_rows,"m.mch_trans"] <- NA # unknown 
flds[tmp_idx_rows,"m.mch_exp"]   <- NA # unknown 
flds[tmp_idx_rows,"m.stage_x"]   <- NA # unknown 
flds[tmp_idx_rows,"m.stage_y"]   <- NA # unknown 
flds[tmp_idx_rows,"m.stage_z"]   <- NA # unknown 

# check again:
tmp <- flds %>% filter(e.folder %in% c("2015_0330_run2_1771wt_15fields_60min_day4_7hrAdeNPG","2015_0821_run1_2149rev7"))

save(flds,file="../rdata_files/flds_with_metadata.RData")


# Re-randomizing to inc Rev7 data on 12/27 --------------------------------

# import the fields + metadata:
load("../rdata_files/flds_with_metadata.RData")
tmp <- flds %>% filter(!complete.cases(.))

# wee need to also import the exps data from the excel file, bc this contains genotype info:
xl <- "../raw_scoring_data/2016_blinded_analysis.xlsx"
exps <- read_excel(xl, sheet = "exps")
exps <- exps %>% 
  # remove rows containing only NA
  filter( rowSums(is.na(.)) < ncol(.) ) %>% 
  # only keep rows with e.analyse_aug16 == TRUE
  mutate(e.analyze_aug16 = as.logical(e.analyze_aug16)) %>% filter(e.analyze_aug16)


# merge them here:
exps.plus.flds <- inner_join(exps, flds, by = c("e.folder")) %>% arrange(f.order_aug16)
# this first merge should NOT change the number of rows
stopifnot( identical(nrow(exps.plus.flds), nrow(flds)) )

# make tier_done, tier_1, tier_2, etc. based on priority
tier_done <- exps.plus.flds %>% filter(f.order_aug16 <= 300)
not_done  <- exps.plus.flds %>% filter(f.order_aug16 > 300)

tier_1  <- not_done %>% filter(e.genotype == "rev7")
tier_2  <- not_done %>% filter(e.genotype %in% c("rev7exo1","exo1","crb2rqh1","rev7crb2"))
tier_3  <- not_done %>% filter(e.genotype %in% c("wt","rqh1"))
tier_4  <- not_done %>% filter(e.genotype %in% c("crb2","crb2exo1","rev7rqh1","rev3"))
tier_5  <- not_done %>% filter(e.genotype %in% c("pht1","crb22AQ"))

stopifnot(identical(nrow(rbind(tier_done,tier_1,tier_2,tier_3,tier_4,tier_5)),nrow(flds)))


# combine together and randomize rows based on tier priorities
rand_order_2 <- sample.int(nrow(tier_2))
rand_order_3 <- sample.int(nrow(tier_3))
new301_to_400 <- bind_rows(tier_1,slice(tier_2,rand_order_2[1:(100-nrow(tier_1))])) %>% sample_frac(1)
new401_to_500 <- slice(tier_2,rand_order_2[ 1:100 + (100-nrow(tier_1))  ]) %>% sample_frac(1)
new501_to_600 <- bind_rows(slice(tier_2,rand_order_2[ (100-nrow(tier_1)+100+1) : nrow(tier_2)  ]),
                           slice(tier_3,rand_order_3[ 1 : (300-nrow(tier_1)-nrow(tier_2))  ])) %>% sample_frac(1)
new602_to_endoftier3 <- slice(tier_3,rand_order_3[(300-nrow(tier_1)-nrow(tier_2)+1) : nrow(tier_3)]) %>% sample_frac(1)
require(lubridate)
new_flds <- bind_rows(tier_done,
                     new301_to_400,
                     new401_to_500,
                     new501_to_600,
                     new602_to_endoftier3,
                     sample_frac(tier_4,1),
                     sample_frac(tier_5,1)) %>%
            select_(.dots = names(flds)) %>%
            mutate(m.t_int = duration(m.t_int), m.t_tot = duration(m.t_tot))

stopifnot(identical(flds,arrange(new_flds,f.order_aug16))) # problem is the time columns are changing class!!!

# make a new column called f.order_dec16
new_flds <- new_flds %>% mutate(f.order_dec16 = row_number()) %>% select(1:f.order_aug16,f.order_dec16,everything())

# RData file of flds that will be used in other scripts
flds <- new_flds
save(flds,file="../rdata_files/flds_with_metadata_newDec16.RData")

# .tsv file that will be used by imageJ macros:
flds %>% select(e.folder,f.file_prj,f.rand_lab,f.order_aug16,f.order_dec16) %>% 
         write.table("../ij_macros/flds_blinding_file_dec16.tsv",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)




# Randomizing image analysis ~Aug 2016 ------------------------------------

# Import and clean .xlsx sheet "exps" -------------------------------------
exps <- read_excel("../2016_blinded_analysis.xlsx", sheet = 1)

# remove rows with all NAs and anything that is 1 for col analyse_aug16
exps <- exps %>% filter( rowSums(is.na(.)) < ncol(.) ) %>% 
                 mutate(e.analyse_aug16 = as.logical(e.analyse_aug16)) %>% filter(e.analyse_aug16)

# remake categorical variables as factors (sort the oder of the factors for genotype)
tmp_levels <- c("wt","exo1","rqh1",
                "crb2","crb2exo1","crb2rqh1",
                "rev7","rev7exo1","rev7rqh1","rev7crb2",
                "pht1","rev3","crb22AQ")
exps <- exps %>% mutate(run = factor(run), genotype = factor(genotype,tmp_levels))

# Information about the experiments so far --------------------------------

# how many total fields, imaging replicates, and biological replicates do we have for each??
exps %>% group_by(genotype) %>% 
         summarise(fields = sum(n_fields), 
                   img_reps = length(n_fields), 
                   biol_reps = n_distinct(date)) %>% 
         arrange(desc(fields))

# Make a new df, "flds", with each field as a separate row -------------------------

path.root <- file.path( getwd(),"..", "PRJs_only" );

get_dvs <- function(thisfolder) {
  filesall <- dir(file.path(path.root,thisfolder))
  filesprjs <- str_subset(filesall,"^[^~\\.].*\\.dv$")
}

# each row in this df will be the .dv file from a single filed, plus all columns
# from exps as well!
flds <- exps %>% group_by(folder) %>% 
                 do(data.frame(prj_file = get_dvs(.$folder))) %>% 
                 full_join(exps,.,by="folder")

# search PRJ file names for a number "..._##_...dv" as close to the end as possible
# str_match(...)[,2] returns the 1st capturing group of the regexp (## only!)
flds <- flds %>% mutate( field = as.numeric(  str_match(prj_file, "_(\\d{1,2})_(?!.*_)")[,2]  ))

# how many fields are we missing PRJ files for??
flds %>% group_by(folder,n_fields) %>% summarize(n_files=n(),diff=first(n_fields)-n()) %>% print(n=nrow(.))

# check how many fields we have for each genotype now that we know which files are missing
flds %>% group_by(genotype) %>% 
         summarise(fields = n_distinct(prj_file), 
                   img_reps = n_distinct(folder), 
                   biol_reps = n_distinct(date)) %>% 
         arrange(desc(fields))

# generate unique strings to anonymously label each PRJ file --------------
rand_lab <- character(0)
while (length(rand_lab) < nrow(flds)) {
  rand_lab <- stri_rand_strings(nrow(flds) + 100, 4, pattern="[a-z0-9]")
  rand_lab <- unique(rand_lab)
}
flds$rand_lab <- rand_lab[1:nrow(flds)]

# copy to clipboard (past into excel sheet ... this mangles some of the random char strings!!)
flds %>% select(folder,prj_file,rand_lab) %>% write.table("flds_temp.tsv",sep="\t",row.names=FALSE,col.names=TRUE)

# save to temporary .xlsx ... then copy paste into my workbook!
# I couldn't get it to work this way, so I went ahead and just did the above command instead!
require(xlsx)
write.xlsx2(flds,"fldx_temp.xlsx",row.names=FALSE,showNA=FALSE)






