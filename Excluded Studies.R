### Excluded Studies ###
library("readxl")

site.info <- read_excel(path = "data/160210_Altadapt_Study.xlsx", sheet = 1, col_names = TRUE)

# Exclude empty rows
site.info <- site.info[,-c(33, 35:39)]

site.info <- site.info %>%
  as.tibble() %>% 
  filter(!is.na(SEQ1)) %>% 
  mutate(ID = paste(SiteID, species, sep="_")) %>% # create ID for species within StudyID
  arrange(ID) %>% 
  select(StudyID, SiteID, ID, species, studysite, data_type, Coord, coor.var, LATLONG_extr, excluded, why_exclude, COMMENT) 
  

# Count the number of studies that were screened
site.info %>% 
  group_by(StudyID, excluded) %>% 
  summarize(n = n())
# 216 studies

included.studies <- dat2 %>% 
  distinct(StudyID, ID)

transplants <- site.info %>% 
  filter(data_type == "transplant") %>% 
  distinct(StudyID, ID)

# How many studies were screened and not used
not.used <- site.info %>% 
  anti_join(included.studies, by = c("StudyID", "ID"))


not.used %>% 
  # Excluded because only 2 pops, cannot calc slopes
  mutate(why_exclude = ifelse(StudyID %in% c("Bast2015JPlE", "Cast2013PlED", "Fetc2000Biot", "Frei2014PloO", "Gale1991Evol"), "only_2_pops", why_exclude)) %>% 
  # Excluding studies with mean.dist.km < 1km
  mutate(why_exclude = ifelse(StudyID %in% c("Gonz2009JEco", "Haut2009JPlE", "Sund1995Scan", "Fisc2011AlpB", "Gime2007AnBo", "Mari1993aEcoR", "Rice1991Oeco", "Saen2011Agro", "Saen2013APhP", "Sone2007Bryo"), "small_elev", why_exclude)) %>% 
  # Excluding studies with dist.km > 1000km
  mutate(why_exclude = ifelse(StudyID %in% c("Ginw2004aSil", "Ginw2004bSil", "Gomo2011AnPS", "Larw2010Silv", "Nels1967NePh"), "large_lat", why_exclude)) %>% 
  # no variation in coordinates or no data available, only Latitude available
    mutate(why_exclude = ifelse(StudyID %in% c("Rehf1983aCJFR", "Rehf1984CJFR", "Rehf1988SiGe", "Rehf1989ForE", "Rehf1990BoGa", "Gaut1998NewP", "Grav1988NewP", "Gure1992Func", "Gure1992Gene", "Hamr1976Theo", "Emer1994IJPl", "Ishi2013PlBi", "Kard2014RSOS", "Gure1988AmJB", "Hove2003NewP", "Plue2005FunE", "Rehf1979Amer", "Rehf1983bCJFR", "Rehf1986USDA", "Stin2004Amer", "Vick1983Grea", "Sund1995Scan", "Poll2009DiDi", "Rehf1994CaJB", "Adam2007Indu", "Ares2000IJPl", "Enss2015AmJB", "Fait2000PIHS", "Fish2007AJBo", "Haid2011PlED", "Hook1990Biom", "Kim2011AmJB", "Lems1971Ecol", "Link2003GlCB", "McKa2001PRSB", "Melc2000IJPS", "Rehf1986ForS", "Rehf1986USDA", "Rehf1991CJFR", "Rehf1993AmJB", "Rehf1994CJFR", "Rein2011Tree", "Sout2011AnBo", "Stin2005Arct", "Vera1997Plan", "Vita2009Cana", "XieC1995CJFR", "Rosn2003JRMa", "Weng1997Phot"), "no_data", why_exclude)) %>% 
  # no data collected, because traits not used
  filter(!StudyID %in% c("Vita2014FunE", "Saen2013APhP", "Zhan2015PlSB")) %>% 
  # Exclude transplants and studies where some data was used
  filter(!StudyID %in% c("Halb0001subm", "Hall1990CJPl","Lege2009MolE", "Voli2002Biol", "Haid2012Oeco", "Vita2013Oeco")) %>% 
  filter(StudyID != "NA", !why_exclude %in% c("duplicate", "NA")) %>% 
  distinct(why_exclude, StudyID) %>% 
  group_by(why_exclude) %>% 
  summarise(n = n(), percent = n * 100 / 201)
  
  