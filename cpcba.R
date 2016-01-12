library(phyloseq)
library(ape)
library(ggplot2)

#load otu table from mothur shared file (edited manually)
otus <- as.matrix(read.table("~/ssd/mothur/0.03.shared", fill = T, header = T, row.names = 1))

#load metadata
sample.data <- read.table("~/Work/Anderson_Lake/key3.csv", sep = "\t",header = F, row.names = 1, col.names = c("sample","id","lake","date","anatoxin"))
sample.data$anatoxin <- as.factor(sample.data$anatoxin)
sample.data$date <- as.Date(sample.data$date, format = "%m/%d/%Y")
#remove samples in metadata that are not in the otu table
sample.data <- sample.data[rownames(otus),]

#make phyloseq object
cpc <- phyloseq(otu_table(otus, taxa_are_rows = F), sample_data(sample.data))

#remove otus with fewer than 20 sequences in total
cpc <- prune_taxa(taxa_sums(cpc)>=20,cpc)

#write the pruned otus to a file to make a tree
write.table(taxa_names(cpc),file = "~/ssd/mothur/over20sum.otus",quote = F,row.names = F, col.names = F)

#load a tree made from FastTree using the aligned fastas from the pruned otus (using the get.otureps function in mothur)
tree <- read.tree("~/ssd/mothur/over20sum.tre")

#add the tree to the phyloseq object
cpc <- merge_phyloseq(cpc, tree)

#remove the control samples from the phylseq object and create a new one
cpc2 <- prune_samples(!sample_names(cpc)%in%c("H12","H11"),cpc)

#remove otus and samples with fewer than 1000 counts in any one sample
cpc3 <- prune_taxa(colnames(get_taxa(cpc))[which(get_taxa(cpc)>=1000,arr.ind = T)[,2]],cpc)
cpc3 <- prune_samples(rownames(get_taxa(cpc3))[which(get_taxa(cpc3)>=1000, arr.ind = T)[,1]],cpc3)

#set default ggplot2 theme
theme_set(theme_bw())

#plot heatmap of lakes ordered by weighted unifrac of community composition
plot_heatmap(cpc2, sample.label = "lake")

#plot bar graph of otu abundances in each lake, colored by anatoxin pcr detection level
plot_bar(cpc3, x = "lake", fill = "anatoxin")

#plot tree of otus with anatoxin measurements on each leaf
plot_tree(cpc, label.tips = "taxa_names", color = "anatoxin", size = "abundance", base.spacing = 0.05, sizebase = 2, min.abundance = 1000)

#MDS plot of weighted unifrac distances between lake otu communities
mds <- ordinate(cpc2, "MDS", distance = "wunifrac")
q <- plot_ordination(cpc2, mds, shape = "anatoxin", color = "lake")
q <- q + ggtitle(paste("MDS using distance method ", i, sep = "")) + geom_text(aes(label = id, color = lake, check_overlap = T), size = 3, nudge_y = 0.01)
q
