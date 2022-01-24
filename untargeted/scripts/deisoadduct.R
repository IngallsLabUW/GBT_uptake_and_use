# Remove duplicates (choose biggest) and fill in missing
filled_peaks <- raw_peaks %>%
  arrange(desc(into)) %>%
  group_by(feature, filename) %>%
  slice(1) %>%
  ungroup() %>%
  complete(feature, nesting(filename, sample)) %>%
  group_by(feature) %>%
  mutate(across(starts_with("rt"), ~ifelse(is.na(.x), median(.x, na.rm = TRUE), .x))) %>%
  mutate(mz=ifelse(is.na(mz), mean(mz, na.rm=TRUE), mz))
first_min <- filled_peaks %>%
  group_by(feature) %>%
  summarise(rtmed=median(rt)) %>%
  filter(rtmed<60)
filled_peaks <- filled_peaks %>%
  filter(!feature%in%first_min$feature)

# Find isotopes for all peaks
message("Finding isotopes for all peaks...")
peak_envelopes <- filled_peaks %>%
  split(.$filename) %>% 
  bplapply(FUN = findIsoAdduct, xdata=xdata_filled,
           grabSingleFileData=grabSingleFileData, checkPeakCor=checkPeakCor, 
           pmppm=pmppm, trapz=trapz, polarity=polarity, mzxml_path=mzxml_path) %>%
  rbindlist() %>%
  `rownames<-`(NULL) %>% arrange(feature)

# Identify adducts and isotopes within the peaks
message("Identifying adducts and isotopes...")
is_iso_adduct <- bplapply(split(filled_peaks, filled_peaks$filename),
         FUN = isIsoAdduct, xdata=xdata_filled,
         grabSingleFileData=grabSingleFileData, checkPeakCor=checkPeakCor, 
         pmppm=pmppm, trapz=trapz, polarity=polarity, mzxml_path=mzxml_path) %>%
  do.call(what = rbind) %>% as.data.frame()

peakshapematch <- is_iso_adduct %>%
  group_by(feature) %>%
  summarise(across(contains("match"), median)) %>%
  as.data.frame()

peakareamatch <- lapply(unique(is_iso_adduct$feature), function(i){
  feature_areas <- is_iso_adduct[is_iso_adduct$feature==i,]
  area_cols <- grep(pattern = "area$", names(feature_areas), value = TRUE)[-1]
  sapply(area_cols, function(x){
    suppressWarnings(cor(feature_areas$M_area, feature_areas[[x]]))
  })
}) %>% 
  do.call(what=rbind) %>% 
  `[<-`(is.na(.), 0) %>% 
  as.data.frame(stringsAsFactors=FALSE) %>%
  mutate(feature=unique(is_iso_adduct$feature)) %>%
  select(feature, everything()) %>%
  arrange(feature)



addiso_features <- peakshapematch[
  which(rowSums(peakareamatch[,names(peakareamatch)!="feature"]>
                  area_remove_threshold&
                  peakshapematch[,names(peakshapematch)!="feature"]>
                  shape_remove_threshold)>=1),
]
addiso_features$adduct_type <- addiso_features %>%
  split(seq_len(nrow(.))) %>%
  sapply(function(i){
    names(which.max(i[-1]))
  }, USE.NAMES = FALSE) %>%
  gsub(pattern = "_match", replacement = "")

safe_features <- sapply(not_addisos, function(cmpd_data){
  filled_peaks %>%
    filter(mz%between%pmppm(cmpd_data[["mz"]], 10)) %>%
    filter(rt%between%(cmpd_data[["rt"]]+peak_rt_flex*c(-1, 1))) %>%
    pull(feature) %>%
    unique()
}) %>% unlist()
addiso_features <- addiso_features %>%
  filter(!feature%in%safe_features) %>%
  select(feature, adduct_type)


# Combine for prettiness' sake
filled_peaks <- filled_peaks %>%
  left_join(peak_envelopes %>% select(feature, filename, M_area))
