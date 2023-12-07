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
library(circlize)
library(BioCircos)
library(SynExtend)
#library(ggmsa) 
#BiocManager::install("ggmsa")
#Phylogeny####
# Your list of accessions
accession_list <- c("KF986481","AY580003", "MG320838", "MN184018","MW124877", "JQ306054", "LC053597", 
                    "HQ241545","JQ306102", "MN184090", "MK916102", "MG320510", "MG321135", "AB857346",
                    "AB114222", "EU664339","AB857347", "MK308325", "KP255317", "MN184051", "KU906044",
                    "MH087706", "HM638048","EF683576")
#True crabs vs false crabs matched with accession list
true_false <- c("Root", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE
                , TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE)
#get sequences and write file
if(file.exists("crabs.fasta")== FALSE){
  for (i in accession_list){
test <- entrez_fetch(db = "nucleotide",id = i, rettype = "fasta")
test <- test %>% str_replace("\n","~") %>% str_remove_all("\n") %>% str_split("~") %>% unlist()
write_lines(test, "crabs.fasta", append = TRUE)
}
}
#read back in file
fasta <- readFasta("crabs.fasta")
# make dataframe with names and sequences
id <- fasta@id %>% as.character %>% unlist
df <- data.frame(accession = id %>% str_split(" ") %>% map_chr(1),
                 species = paste(id %>% str_split(" ") %>% map_chr(2),id %>% str_split(" ") %>% map_chr(3), sep = " "),
                 true = true_false,
                 sequence = fasta@sread)
#align the sequences and add the names and true vs fake crab
align <- AlignSeqs(fasta@sread)

names(align) <- df$species
#make and root tree
dist_mat <- dist.dna(as.DNAbin(align))
tree <- nj(dist_mat)
rooted <- ape::root(tree, outgroup = 1, resolve.root = TRUE)
rooted$edge.length <- rep(.1, Nedge(rooted)) 
#Export tree to a file (e.g., Newick format)
#write.tree(rooted, file = "phylogenetic_tree.nwk")

p$data$isTip
#plot tree
p <- ggtree(rooted,layout = "rectangular") + 
  theme(plot.margin = margin(0,0,0,0))+
  coord_cartesian(xlim = c(0,2.1))
p2 <- p + 
  geom_tiplab(align = TRUE, aes(color = df$true), data= p$data %>%dplyr::filter(isTip ==TRUE))+
  ggtitle("Phylogeny of Crabs")+
  theme(legend.title = "Key")

p2
ggsave("crab_phylogeny.png")
#Synteny####
test_string_set <- fasta2AAStringSetlist("../../../crab_full_genome/")
seqs <- c(Birgus_latro = "../../../crab_full_genome/birgus_latro/ncbi_dataset/data/GCA_018397915.1/filtered.fna",
Chionecetes_opilio = "../../../crab_full_genome/chionoecetes_opilio/ncbi_dataset/data/GCA_016584305.1/filtered.fna"
)

db <- "../../../crab_full_genome/new_db.fasta"

for (i in seq_along(seqs)) {
  Seqs2DB(seqs[i], "FASTA", db, names(seqs[i]))
}

syn <- FindSynteny(db,minScore=25)
plot(syn, horizontal=TRUE,width=.6, colorRamp = colorRampPalette(c("red", "orange",
                                                                               "yellow", "green",
                                                                               "blue", "purple")))
class(syn)
ggplot(syn)
print(syn)
pairs(syn)

chordDiagram(syn)
class(syn)
data <- data.frame(syn[[2]])

chordDiagram(data)
