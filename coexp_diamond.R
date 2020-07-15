#!/usr/bin/env Rscript

## Discover disease modules from the RNA-seq data

#Load packages
suppressPackageStartupMessages(library(RankProd))
suppressPackageStartupMessages(library(gmp))
suppressPackageStartupMessages(library(CoExpNets))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(filesstrings))

#Define functions

seed_modules = function(module_color, modules = coexp_net$moduleColors){
  module_genes = names(coexp_net$moduleColors)[coexp_net$moduleColors == module_color]
  return(module_genes)
}


network_mapper = function(genes, nw = ppi){
  merge_interactorA = merge(nw, genes, by.x = 1, by.y = 1)
  merge_interactorB = merge(merge_interactorA, genes, by.x = 2, by.y = 1)
  return(merge_interactorB)
}

#Set input parameters
option_list <- list(
  make_option(c("-i","--input"), type="character", default=NULL,
              help="path to the RNAseq expression matrix txt file, \t [default %default]",
              dest="input_filepath"),
  make_option(c("-w","--network"), type="character", default=NULL,
              help="path to the PPI network, \t [default %default]",
              dest="network_filepath"),
  make_option(c("-s", "--seeds"), type="character", default=NULL,
              help="path to the seed genes for the network construction [default %default]",
              dest="input_s_filepath"),
  make_option(c("-o","--output"), type="character", default="coexp_diamond",
              help="output filename prefix [default %default]",
              dest="output_filename"),
  make_option(c("-a", "--alpha"), type="numeric", default=1,
              help="weight for the seed genes [default %default]",
              dest="a"),
  make_option(c("-n", "--number_interactors"), type="numeric", default=300,
              help="Number of interations performed by DIAMOnD, equal to number of interactors in the disease network [default %default]",
              dest="n"))

parser <- OptionParser(usage = "%prog -i exp_matrix.txt -s seeds.txt -w ppi_network.txt [options]",option_list=option_list)
param = parse_args(parser)

cat(paste0("PARAMETERS:\ninput: ", param$input_filepath,"\n"))
cat(paste0("network: ", param$network_filepath,"\n"))
cat(paste0("output: ", param$output_filename,"\n"))
cat(paste0("seed genes: ", param$input_s_filepath,"\n"))
cat(paste0("alpha: ", param$a,"\n"))
cat(paste0("Number of DIAMOnD iterations: ", param$n,"\n"))

#Read input files
edata = read_delim(param$input_filepath, "\t", escape_double = FALSE, trim_ws = TRUE)
genes = as.vector(unlist(edata[, 1]))
edata = edata[, -1]
rownames(edata) = genes
edata = as.data.frame(t(edata))

ppi = read_delim(param$network_filepath, "\t", escape_double = FALSE, trim_ws = TRUE)

seeds = read_delim(param$input_s_filepath, "\t", escape_double = FALSE, trim_ws = TRUE)
seeds = as.vector(unlist(seeds[, 1]))

a = param$a
n = param$n

o = param$output_filename

#Create output directories

dir.create(paste0("./", o, "/input_diamond/"), recursive = TRUE)
dir.create(paste0("./", o, "/output_diamond/"), recursive = TRUE)

#Build co-expression networks

print("Build co-expression networks")
coexp_net = CoExpNets::getDownstreamNetwork(tissue="tissue",
                                         n.iterations=20,
                                         net.type = "signed",
                                         debug=F,
                                         expr.data=edata,
                                         job.path=paste(o, "/CoExpNets/"))


nmodules = dim(coexp_net$MEs)[2]
smodules = coexp_net$moduleColors[which(names(coexp_net$moduleColors) %in% seeds)]

#Print results of module selection
print(paste0(nmodules, " modules have been identified. Seed genes are present in ", length(unique(smodules)), " modules."))

module_genes = sapply(unique(smodules), seed_modules)
names(module_genes) = unique(smodules)

module_networks = lapply(module_genes, network_mapper)
diamond_input = paste0("./",o, "/diamond_input/")

for (i in 1:length(module_networks)){
  myfile <- paste0( o, "/input_diamond/", names(module_networks)[i], ".txt")
  write.table(module_networks[[i]], myfile , 
              quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
}


print("Run DIAMOnD")
o_diamond = paste0(o, "/output_diamond/")
i_diamond = paste0(o, "/input_diamond/")

for (i in 1:length(module_networks)){
  system(paste('python2 DIAMOnD.py', paste0(i_diamond, names(module_networks)[i], ".txt"), param$input_s_filepath,
               a, n, paste0(o_diamond, names(module_networks)[i], "_",n,"_nodes.txt"), sep = " "), wait = FALSE)
}


