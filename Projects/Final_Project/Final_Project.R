library(tidyverse)
library(ape)
library(phangorn)
library(janitor)
library(rentrez)
library(Biostrings)
library(seqinr)
library(ShortRead)
library(DECIPHER)
library(ggtree)
# Your list of accessions
accession_list <- c("KF986481","AY580003", "MG320838", "MN184018","MW124877", "JQ306054", "LC053597", 
                    "HQ241545","JQ306102", "MN184090", "MK916102", "MG320510", "MG321135", "AB857346",
                    "AB114222", "EU664339","AB857347", "MK308325", "KP255317", "MN184051", "KU906044",
                    "MH087706", "HM638048","EF683576")
true_false <- c("Root", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE
                , TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE)
i <- "KF986481"
if(file.exists("crabs.fasta")== FALSE){
  for (i in accession_list){
test <- entrez_fetch(db = "nucleotide",id = i, rettype = "fasta")
test
test <- test %>% str_replace("\n","~") %>% str_remove_all("\n") %>% str_split("~") %>% unlist()
test
write_lines(test, "crabs.fasta", append = TRUE)
}
}
fasta <- readFasta("crabs.fasta")
id <- fasta@id %>% as.character %>% unlist
base::paste(as.character(id[[1]][2:3]), sep = " ")

df <- data.frame(accession = id %>% str_split(" ") %>% map_chr(1),
                 species = paste(id %>% str_split(" ") %>% map_chr(2),id %>% str_split(" ") %>% map_chr(3), sep = " "),
                 true = true_false)


idk <- AlignSeqs(check)
#aligned_maybe <- as(idk, "DNAStringSet")
idk
dist_mat <- dist.dna(as.DNAbin(idk))
tree <- nj(dist_mat)

plot(tree)  # Basic tree plot

# Export tree to a file (e.g., Newick format)
write.tree(tree, file = "phylogenetic_tree.nwk")

ggtree(tree,layout = "slanted") + geom_tiplab()+
  geom_rootedge(rootedge = )


