require(readxl)
require(stringi)
require(stringr)
require(lubridate)
require(magrittr)
require(tidyr)
require(dplyr)


#### Still got to FULL get the metadata for 191 383 471 693 764?
#### try to get it from original .dv files (not PRJ)


# Hardcoded section: ------------------------------------------------------

xl <- "../raw_scoring_data/2016_blinded_analysis.xlsx"

genotype.factor.order <- c("wt","exo1","rqh1",
                           "crb2","crb2exo1","crb2rqh1",
                           "rev7","rev7exo1","rev7rqh1","rev7crb2",
                           "pht1","rev3","crb22AQ")


# Import from 2016_blinded_analysis.xlsx & flds_with_metadata.RData --------
# do a few minor tweaks and validations that need to happen before
# merging cells with exps and flds

# experiments:
exps <- read_excel(xl, sheet = "exps")
exps <- exps %>% 
  # remove rows containing only NA
  filter( rowSums(is.na(.)) < ncol(.) ) %>% 
  # only keep rows with e.analyse_aug16 == TRUE
  mutate(e.analyze_aug16 = as.logical(e.analyze_aug16)) %>% filter(e.analyze_aug16)
  # cannot have any missing data in these columns:
  exps.key.cols <- select(exps, e.genotype,e.strain,e.date,e.run,e.folder,e.n_fields,e.t_start)
  stopifnot( all(!is.na(exps.key.cols)) )
  # must all be unique combinations for these key columns
  stopifnot( nrow(exps) == nrow(distinct(exps.key.cols)) )
  # must ONLY use genotypes we expect; all genotypes must be represented!
  stopifnot( setequal(exps$e.genotype, genotype.factor.order) )


# fields:
# this has already been created and validated by "parse_dv_log_files_for_metadata.R" and "randomizing_script.R"
load("../rdata_files/flds_with_metadata_newDec16.RData")
  

# cells:
# (force all cols to be text bc of "?" and "Inf" ... convert to numeric later)
cells <- read_excel(xl, sheet = "cells", col_types = c("numeric","text","numeric",rep("text",8)))

# cannot have empty rows in cells
stopifnot( all(rowSums(is.na(cells)) < ncol(cells)) )

# add a new column 
cells <- cells %>% 
  # add a new column for the row number in the excel sheet
  mutate(c.xlsx_row = row_number() + 1) %>% 
  select(starts_with("f"),c.xlsx_row,everything())

# validate
cells.validate <- cells %>% 
                  filter(!is.na(f.order_dec16)) %>%
                  mutate(v.is.consec = (f.order_dec16 - lag(f.order_dec16,default=0)) == 1,
                         v.all.na    = select(.,starts_with("c."),-contains("xlsx_row")) %>% {ncol(.)},
                         v.is.valid  = v.all.na | (c.num_line == 1 & !is.na(c.num_line)))
if (any(!cells.validate$v.is.consec)) {
  stop("The following rows of the 'cells' tab may not be consecutive\n\t",
       paste(cells.validate %>% filter(!v.is.consec) %>% .$c.xlsx_row,collapse=","))
}
if (any(!cells.validate$v.is.valid)) {
  stop("The following rows of the 'cells' tab have problems\n",
       "the f.order_ number might be in the wrong row for the following excel rows:\n\t",
       paste(cells.validate %>% filter(!v.is.valid) %>% .$c.xlsx_row,collapse=","))
}


# now that we've validated everything, fill in missing field numbers in 1st col
cells <- cells %>% fill(f.order_dec16)



# cell length measurments
cell_len <- read_excel(xl, sheet = "cell_len")
cell_len <- 
  cell_len %>%
  mutate(c.len_rad52_on = as.numeric(c.len_rad52_on),
         c.len_lacO_off = as.numeric(c.len_lacO_off))


# Manually go through and make c.notes consistent here --------------------


# use this command to help manually go through f.notes
# used lots of find/replace in the .xlsx to clean these up
tmp <- cells %>% filter(!is.na(f.notes)) %>% .$f.notes %>% 
                  str_split("; ") %>% 
                  unlist() %>% tibble(f.notes=.) %>%
                  group_by(f.notes) %>% mutate(n = n()) %>% filter(row_number() == 1) %>%
                  arrange(desc(n),f.notes)

# use this command to help manually go through c.notes
# used lots of find/replace in the .xlsx to clean these up
tmp <- cells %>% filter(!is.na(c.notes)) %>% .$c.notes %>% 
  str_split("; ") %>% 
  unlist() %>% tibble(c.notes=.) %>%
  group_by(c.notes) %>% mutate(n = n()) %>% filter(row_number() == 1) %>%
  arrange(desc(n),c.notes)



# Merge cells, exps, flds, and lengths ---------------------------------------------

exps.plus.flds <- inner_join(flds, exps, by = c("e.folder"))
# this first merge should NOT change the number of rows
stopifnot( identical(nrow(exps.plus.flds), nrow(flds)) )
# also should preserve the order of f.order_dec17
stopifnot(
  exps.plus.flds$f.order_dec16 %>% unique() %>% diff() %>% {.==1} %>% all()
)


# create df, a full data frame that includes ALL of the cell data
# but also includes all of the exp/field data that has not yet
# been analyzed (<NA> for all of the cells columns) ...
# that is why this uses full_join!!
df <- full_join(exps.plus.flds,cells,by=c("f.order_dec16"))
# df should have N rows, where 
#   N = numb.of.fields + numb.of.cells.analyzed - numb.of.fields.completed.in.analysis
stopifnot( nrow(df) == nrow(flds) + nrow(cells) - max(cells$f.order_dec16) )
stopifnot(
  df$f.order_dec16 %>% unique() %>% diff() %>% {.==1} %>% all()
)

# add on the c.len_... columns from cell_len to the final df
df <- full_join(df,
          cell_len %>% select(c.xlsx_row,c.len_rad52_on,c.len_lacO_off),
          by = "c.xlsx_row")




# Clean up and validate df, now that it has been fully merged -------------

# validate f.notes and f.last_good_frame
df.f.validate <- df %>% group_by(f.order_dec16) %>% 
                        mutate(v.misplaced = (row_number() != 1) & (!is.na(f.notes) | !is.na(f.last_good_frame)),
                               v.last_frame_bad = !is.na(f.last_good_frame) & (f.last_good_frame > m.t_num))
if (any(df.f.validate$v.misplaced)) {
  stop("The following rows of the 'cells' tab have misplaced info entered\n",
       "in either the f.notes or the f.last_good_frame columns:\n\t",
       paste(df.f.validate %>% filter(v.misplaced) %>% .$c.xlsx_row,collapse=","))
}
if (any(df.f.validate$v.last_frame_bad)) {
  stop("The following rows of the 'cells' tab have too high a value for f.last_good_frame\n\t",
       paste(df.f.validate %>% filter(v.last_frame_bad) %>% .$c.xlsx_row,collapse=","))
}
# ensure that these two columns are the correct class (didn't misplace any charachter in the f.last_good_frame col!)
stopifnot(class(df$f.notes) == "character")
stopifnot(class(df$f.last_good_frame) == "numeric")


# do some easy clean-ups on some columns on the left side of df
# mostly related to the experiment ("e.___") and field ("f.___")
df <- df %>%
  # turn e.genotype column into factor
  mutate(e.genotype = factor(e.genotype,levels=genotype.factor.order)) %>%
  # convert to minutes (Duration)
  mutate(e.t_slide = dminutes(e.t_slide), e.t_start = dminutes(e.t_start) ) %>%
  # fill in missing f.notes and f.last_good_frame
  group_by(f.order_dec16) %>% fill(f.notes) %>% fill(f.last_good_frame) %>% ungroup() %>%
  # if f.last_good_frame is empty (NA) then use the end of the movie as the last good frame
  mutate(f.last_good_frame = if_else(is.na(f.last_good_frame),m.t_num,f.last_good_frame)) %>%
  # make a new column saying whether or not a field has been analyzed yet
  mutate(f.is_analyzed = f.order_dec16 %in% cells$f.order_dec16) %>%
  # make a new column for fields that have been analyzed, but are not usable
  # these are the fields that have only one row in the 'cells' sheet and
  # nothing (NA) in the c.num_line column 
  group_by(f.order_dec16) %>% mutate(f.is_usable = !all(is.na(c.num_line))) %>% ungroup()





# start doing cell-specific ("c.___") validation and cleaning here!


# c.num_line must contain only numbers or "a","b","aa","ab",...  
# for the rows that have been analyzed and where the field was useable
df.c.validate <- df %>% filter(f.is_analyzed) %>% transmute(c.xlsx_row, bad_rows = !str_detect(c.num_line,"^[\\dab]+$") & f.is_usable)
if (any(df.c.validate$bad_rows)) {
  stop("The following rows of the 'cells' tab had problems with the format\n",
       "of the column 'c.num_line' (must be a number or 'a','b','ab',...):\n\t",
       paste(df.c.validate %>% filter(bad_rows) %>% .$c.xlsx_row,collapse=","))
}

# split c.num_line into a column of numbers (c.numb) and lineage (c.line)
# ex:  c.num_line:  NA  1   a   b   aa  ab  2   NA  3   a   b ...
#          c.numb:  NA  1   1   1   1   1   2   NA  3   3   3 ...
#          c.line:  NA  p   a   b   aa  ab  p   NA  p   a   b ...
df <- df %>%  mutate(c.numb = suppressWarnings( as.numeric(c.num_line) ),
                     c.line = if_else(is.na(c.numb), c.num_line, "p")) %>%
              group_by(f.order_dec16) %>% fill(c.numb) %>% ungroup()

# validate that we didn't skip any cell numbers!
df.c.validate <- df %>% filter(f.is_usable) %>% group_by(f.order_dec16) %>% 
                        filter(!duplicated(c.numb)) %>%
                        mutate(v.is.consec = (c.numb - lag(c.numb,default=0)) == 1)
if (any(!df.c.validate$v.is.consec)) {
  stop("The following rows of the 'cells' tab may not have consecutive cell numbers:\n\t",
       paste(df.c.validate %>% filter(!v.is.consec) %>% .$c.xlsx_row,collapse=","))
}


# Only valid strings are allowed in the columns related to fames where events happen:
# c.div_nuc	c.div_cell	c.rad52_on	c.lacO_off	c.rad52_off	c.lacO_on
df.c.validate <- df %>% select(c.xlsx_row,c.div_nuc:c.lacO_on) %>% 
                        mutate_at( vars(c.div_nuc:c.lacO_on), funs( is.na(.) | str_detect(.,"\\?|\\d|Inf") ) ) 
if (  any(rowSums(df.c.validate[,-1]) != ncol(df.c.validate)-1)  ) {
  stop("The following rows of the 'cells' tab had problems with the format\n",
       "in one of the div or rad52 or lacO columns (must be '?', 'Inf', or a number):\n\t",
       paste(df.c.validate %>% filter(rowSums(.[,-1]) != ncol(.)-1)  %>% .$c.xlsx_row, collapse=","))
}

# Now convert the strings in these rows to <dbl>
# !is_usable        --> NA
# is_usable && "?"  --> NA
# is_usable && "0"  --> -Inf
# is_usable && "Inf"--> Inf
# is_usable && NA   --> Inf
# is_usable && "#"  --> as.numeric("#")
div_str_to_dbl <- function(div_col,is_usable) {
  out <- suppressWarnings( as.numeric(div_col) ) # silently coerce "?" to NA
  out[is_usable & is.na(div_col)] <- Inf
  out[out == 0] <- -Inf
  return(out)
}


df <- df %>% 
  # fill in "0"s for nuclear divisions that happened before the 1st frame
  mutate(c.div_nuc = if_else(is.na(c.div_nuc) & !is.na(c.div_cell),"0",c.div_nuc)) %>% 
  # fill in "0"s for rad52 loading events that happened before the 1st frame
  mutate(c.rad52_on = if_else(is.na(c.rad52_on) & !is.na(c.lacO_off),"0",c.rad52_on)) %>% 
  # use the function above to convert from string to double (w/ NA and Inf and -Inf)
  mutate_at( vars(c.div_nuc:c.lacO_on), funs(div_str_to_dbl(.)), is_usable = .$f.is_usable)


# loop through each row and find out the idx of that cell's parent
# and when that cell and nuc were born.
df <- df %>% mutate(c.born_nuc = as.numeric(NA), 
                    c.born_cell = as.numeric(NA),
                    c.parent_xlsx_row = as.numeric(NA),
                    c.parent_has_dsb = FALSE)

for (r in 1:nrow(df)) {
  sis <- df$c.line[r]
  if (is.na(sis)) {
    # leave c.born_nuc, c.born_cell, and c.parent_xlsx_row as just NA
    # leave c.parent_has_dsb as FALSE
    next 
    
  } else if (identical(sis,"p")) {
    # nuc and cell are born before movie starts (-Inf)
    df$c.born_nuc[r] <- -Inf
    df$c.born_cell[r] <- -Inf
    # c.parent_xlsx_row <- NA
    # c.parent_has_dsb <-  FALSE
    next
    
  } else if (df$c.line[r] %in% c("a","b")) {
    tmp_parent <- "p"
    
  } else { # "ab","aa","ba",etc
    # tmp_parent of "aab" = "aa"; of "ba" = "a"
    tmp_parent <- str_sub( df$c.line[r], end=-2)
  }
  
  # search for the row that has this cell's parent
  tmp_parent_row_idx <- max(which(  df$c.line[1:r] == tmp_parent ))
  
  # save this parent's unique c.xlsx_row here, in case I ever
  # need to find this cell's parent again
  df$c.parent_xlsx_row[r] <- df$c.xlsx_row[tmp_parent_row_idx]
  
  # parent's div times, use those as this rows "born" times:
  df$c.born_nuc[r]  <- df$c.div_nuc[tmp_parent_row_idx]
  df$c.born_cell[r] <- df$c.div_cell[tmp_parent_row_idx]
  
  # was there already a DSB in this cell's parent (finite, NA, or -Inf)?
  # was there already a DSB in this cell's parent's parent? (df$c.parent_has_dsb[tmp_parent_row_idx] == TRUE)
  # if so, it should not be counted in the % rad52 totals
  if (df$c.parent_has_dsb[tmp_parent_row_idx] ||
      is.finite(df$c.rad52_on[tmp_parent_row_idx]) || 
      is.na(df$c.rad52_on[tmp_parent_row_idx]) ||
      (df$c.rad52_on[tmp_parent_row_idx] == -Inf)) {
    df$c.parent_has_dsb[r] <- TRUE
  }
  
  # TODO: link sisters?
}

# validate that there are no negative intervals. Ex: c.div_cell-c.div_nuc > 0,  
# c.lacO_off-c.rad52_on > 0, etc ...
# if any of these is TRUE, then it's a problem. NA is ok
df.c.validate <- df %>% filter(f.is_usable) %>%
                        mutate(v.nuc_to_cell0  = (c.born_cell-c.born_nuc) < 0, 
                               v.nuc_to_cell1  = (c.div_cell-c.div_nuc)   < 0,
                               v.cell_to_cell  = (c.div_cell-c.born_cell) < 0,
                               v.nuc_to_nuc    = (c.div_nuc-c.born_nuc)   < 0,
                               v.rad_to_lacoff = (c.lacO_off-c.rad52_on)  < 0,
                               v.rad_to_radoff = (c.rad52_off-c.rad52_on) < 0,
                               v.lac_to_lacon  = (c.lacO_on-c.lacO_off)   < 0,
                               v.nuc_to_radon  = (c.rad52_on-c.born_nuc)  < 0,
                               v.nuc_to_radoff = (c.rad52_off-c.born_nuc) < 0,
                               v.nuc_to_lacoff = (c.lacO_off-c.born_nuc)  < 0,
                               v.nuc_to_lacon  = (c.lacO_on-c.born_nuc)   < 0)

df.c.validate$v.any_neg_int <- df.c.validate %>% select(starts_with("v.")) %>% apply(1,any,na.rm=T)
if (  any(df.c.validate$v.any_neg_int)  ) {
  stop("The following rows of the 'cells' tab had problems with some negative interval\n",
       "that should not be possible (ex: cell div before nuc div, or lacO_off before rad52_on):\n\t",
       paste(df.c.validate %>% filter(v.any_neg_int) %>% .$c.xlsx_row, collapse=","))
}
# if you need to see more info about the issues, make and View() this df:
tmp <- df.c.validate %>% filter(v.any_neg_int) %>%
                         select(starts_with("v.")) %>% 
                         bind_cols(df.c.validate %>% filter(v.any_neg_int) %>% select(c.xlsx_row,c.parent_has_dsb),.)



#  Manually go back and look at specific cells based on c.notes -----------

# use this command to help manually go through c.notes
# used lots of find/replace in the .xlsx to clean these up
tmp <- df %>% filter(!is.na(c.notes)) %>% .$c.notes %>% 
  str_split("; ") %>% 
  unlist() %>% tibble(c.notes=.) %>%
  group_by(c.notes) %>% mutate(n = n()) %>% filter(row_number() == 1) %>%
  arrange(desc(n),c.notes)

# go back and manually look at all the cells that have one or more of these c.notes:
# "unequal div", "only one daughter gets nuc", "diploid", "failed nuc division spindle collapse?", "spindle problems"
# anything unique c.notes from tmp above
tmp_regex <- c("unequal div",
               "only one daughter gets nuc",
               "diploid",
               "failed nuc division spindle collapse?",
               "spindle problems",
               tmp %>% filter(n == 1) %>% .$c.notes) %>%
  str_c(collapse = "|") %>% 
  str_replace_all( c("\\?" = "\\\\?", "\\(" = "\\\\(", "\\)" = "\\\\)") )

tmp <- df %>% mutate(c.notes.matches = str_detect(c.notes,tmp_regex)) %>% 
              filter(c.notes.matches) %>%
              select(c.xlsx_row,f.order_dec16,c.numb,c.line,c.notes)
# write.csv(tmp,"../progress/manual_unequal_div_second_looks.csv")

tmp <- df %>% filter(str_detect(c.notes,"spindle|nuc")) %>%
              select(e.genotype,c.xlsx_row,f.order_dec16,c.numb,c.line,c.notes)


# Drop unused columns and select() to reorder the columns we want ---------

# use select to reorder the columns
df <- df %>% select(c.xlsx_row,f.order_dec16,f.order_aug16,f.is_analyzed,f.is_usable,f.last_good_frame,e.genotype,
                    c.numb,c.line,
                    c.born_nuc,c.born_cell,c.rad52_on,c.lacO_off,c.div_nuc,c.div_cell,
                    c.rad52_off,c.lacO_on,
                    c.len_rad52_on,c.len_lacO_off,
                    c.parent_xlsx_row, c.parent_has_dsb,
                    c.notes,
                    e.date,e.run,e.folder,f.file_prj,f.rand_lab,f.notes,
                    starts_with("e."),
                    starts_with("m."),
                    everything())

# FINALLY, drop a few cols that we no longer need
df <- df %>% select( -c(e.for_Robert_sept16, e.analyze_aug16, c.num_line) )


# Save the entire data frame here -----------------------------------------

up_to <- df %>% filter(f.is_analyzed) %>% .$f.order_dec16 %>% max() %>% as.character()

save(df,file=str_c("../rdata_files/2016_binded_analysis_upto",up_to,".RData"))



############# WARNING!
# several dplyr commands coerce Durations to numeric (even if you're not acting on the Duration col directly)
# Raised this issue on github here: https://github.com/hadley/dplyr/issues/2214


################ useful for looking at relavent columns:
# spot check: show all lineages that contain one or more DSB events:
df %>% group_by(f.order_dec16,c.line) %>% filter(!all(c.rad52_on == Inf) || !all(c.lacO_off == Inf)) %>% 
  select(e.genotype, starts_with("c"), -c.notes) %>% print(n=Inf)

df %>% filter(f.is_analyzed) %>%
       select(e.genotype,e.t_start,f.order_dec16,f.is_usable,f.notes:c.notes) %>% print(n=100)

z <- df %>% filter(f.is_analyzed) %>%
       select(f.order_dec16, e.genotype, f.notes, f.last_good_frame, m.t_num, 
              c.numb, c.line, c.parent_xlsx_row, c.parent_has_dsb,
              c.born_nuc, c.born_cell, c.div_nuc, c.div_cell, 
              c.rad52_on, c.lacO_off, c.rad52_off, c.lacO_on, c.notes)

# TODO!!!
# ALL cell lineages where there's a DSB, then a division ... maybe clean these up somehow to 
# differentiate between div with one cell with Rad52 only vs both with Rad52???
# anyway, I would need to go through manually to do this!
z <- df %>% filter(f.is_usable) %>% group_by(f.order_dec16,c.numb) %>%
            filter(any(!(c.rad52_on %in% Inf) & !(c.div_nuc %in% Inf))) %>%
            select(c.xlsx_row,
                   f.order_dec16, e.genotype, f.notes, f.last_good_frame, m.t_num, 
                   c.numb, c.line, c.parent_xlsx_row, c.parent_has_dsb,
                   c.born_nuc, c.born_cell, c.div_nuc, c.div_cell, 
                   c.rad52_on, c.lacO_off, c.rad52_off, c.lacO_on, c.notes) %>%
            ungroup()
      
