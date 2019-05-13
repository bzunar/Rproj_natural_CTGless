library(biomartr)
library(Biostrings)
library(coRdon)
library(magrittr)
library(tidyverse)

# parse Scheffersomyces stpitis .keg file (from KEGG)
kegg <- readr::read_lines(file.path("D:", "RStudio", "ldhCU", "pic00001.keg"))
ss_ko <- "pic00001.keg" %>%
    file.path("D:", "RStudio", "ldhCU", .) %>%
    readr::read_lines() %>%
    {.[startsWith(., prefix = "D")]} %>%
    stringr::str_replace("^D      ", "") %>% 
    as.data.frame() %>%
    tidyr::separate(".", into = c("A", "B"), sep = "\t") %>% 
    tidyr::separate("A", into = c("C", "D"), sep = "; ") %>% 
    tidyr::separate("C", into = c("E", "F"), sep = " ") %>% 
    tidyr::separate("B", into = c("G", "H"), sep = "; ") %>% 
    tidyr::separate("G", into = c("I", "J"), sep = " ") %>% 
    dplyr::select(E, I) %>% 
    dplyr::rename(gene = E, ko = I)

# extract only higly expressed ribosomal genes
ss_ribo <- ss_ko[ss_ko$ko %in% coRdon:::RPKOs, ][, 1]
ss_ribo_ko <- ss_ko[ss_ko$ko %in% coRdon:::RPKOs, ][, 2]

# get Scheffersomyces stipitis genome
is.genome.available(db = "refseq", organism = "Scheffersomyces stipitis")

SS.cds.refseq <- getCDS( db       = "refseq", 
                         organism = "Scheffersomyces stipitis",
                         path     = file.path("_ncbi_downloads","CDS"))

Ss_CDS <- read_cds(file     = SS.cds.refseq, 
                   obj.type = "Biostrings")

# rename the genes
names(Ss_CDS) %<>%
    str_match("\\[locus_tag=(PICST_[0-9]+)\\]") %>%
    .[, 2]

# write highly expressed ribosomal genes to file
writeXStringSet(Ss_CDS[ss_ribo], 
                file = file.path("D:", "RStudio", "ldhCU", "ss_ribo.fasta"))

################################################################################
### codon part

# which genes do not contain CTG codon
zeroCTG <- readSet(file = file.path("D:", "RStudio", "ldhCU", 
                                    "seqdump-ldhD-refseq_representative_genomes.fasta")) %>% 
    codonTable() %>%
    .@counts %>% 
    .[, "CTG"] == 0

# calculate MELP
a <- readSet(file = file.path("D:", "RStudio", "ldhCU", 
                              "seqdump-ldhD-refseq_representative_genomes.fasta")) %>%
    .[zeroCTG] %>%
    codonTable()
b <- readSet(file = file.path("D:", "RStudio", "ldhCU", "ss_ribo.fasta")) %>%
    codonTable()
b <- setKO(b, ss_ribo_ko)
kos <- c(getKO(a), getKO(b))
rm(a, b)
a <- readSet(file = file.path("D:", "RStudio", "ldhCU", 
                              "seqdump-ldhD-refseq_representative_genomes.fasta")) %>%
    .[zeroCTG]
b <- readSet(file = file.path("D:", "RStudio", "ldhCU", "ss_ribo.fasta"))

d <- append(a, b) %>%
    codonTable()
d <- setKO(d, kos)

d %>% MELP(ribosomal = TRUE) %>% .[1:length(a), ] %>% summary()
d %>% MELP(ribosomal = TRUE) %>% .[1:length(a), ] %>% which.max()
d %>% MELP(ribosomal = TRUE) %>% .[1:length(a), ] %>% order(decreasing = TRUE)
d %>% MELP(ribosomal = TRUE) %>% .[1:length(a), ] %>% sort(decreasing = TRUE)
d[97]
################################################################################
d[d %>% MELP(ribosomal = TRUE) %>% .[1:length(a), ] %>% order(decreasing = TRUE)] %>%
    .[2]




d[1:10] %>% View()

a %>% View()

width(a) %>% summary()
width(a)[97]


a[zeroCTG] %>% View()
a









