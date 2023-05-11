#!/usr/bin/env Rscript
library(optparse)

option_list <- list( 
  make_option(c("-t", "--tree"), action="store", type="character",default=NULL,
              help="Load phylogenetic tree rescaled in time (nexus format)"),
  make_option(c("-a", "--ancestral"), action="store", type="character",default=NULL,
              help="Load ancestral state reconstruction .json file"),
  make_option(c("-d", "--dpi"), action="store",type="integer", default=NULL, 
              help="Number of days post infection at necropsy",
              metavar="number"),
  make_option(c("-w", "--window"), action="store",type="integer", default=30, 
              help="Length of time windows in days [default %default]"),
  make_option(c("-o", "--output"), action="store",type="character", default="Coevol_output", 
              help="Name for the output file")
)

parser=parse_args(OptionParser(option_list=option_list))

library(treeio)
library(tidyverse)
library(tidytree)
library(plyr)
library(dplyr)
library(tidyr)
library(jsonlite)
library(parallel)
library(ape)
library(data.table)
library(Rgraphviz)

numCores=detectCores()
`%notin%` = Negate(`%in%`)
options(warn=-1)

b = read.beast(parser$tree)
b.tree = as_tibble(b)

tree_data = select(b.tree, node, parent, height, label, states) 
tree_data$height = as.numeric(tree_data$height)
tree_data = tree_data %>% mutate(end_dpi = parser$dpi - height)
colnames(tree_data)[5] = 'end_tissue'

# Extract mutations along branches from json file output by the following command in nextstrain:
mut_info = read_json(parser$ancestral)$nodes

# Each branch item in the list contains both the ancestral sequence (at which end??) and the mutations - only want mutations:
mut_info = lapply(mut_info, "[[", "muts")
# Since several mutations can exist at a single branch, this information is in the form of a list, but it is easier
# to deal with dataframes:
mut_info = lapply(mut_info, function(x) {
  x=do.call(rbind, x)
  x=data.frame(muts = x)
  return(x)
})

# Determine frequency of mutations at individual sites across all branches int he tree and filter out sites that exhibit mutations
# across more than one branch
mut_info = data.table::rbindlist(mut_info, idcol="label") %>%
  mutate(site = as.numeric(gsub("[A-Z-](\\d+)[A-Z-]", "\\1", muts))) %>%
  group_by(site) %>%
  mutate(num_branches = length(unique(label))) %>%
  filter(num_branches>=2) %>%
  arrange(label, site)

#obtaining origin dpi and origin tissue
tree_data$origin_dpi=0
tree_data$origin_tissue= ''

for (i in 1:length(tree_data$parent)){
  tree_data$origin_dpi[i]= tree_data$end_dpi[tree_data$parent[i]]
  tree_data$origin_tissue[i]= tree_data$end_tissue[tree_data$parent[i]]
} 

friedrice = merge(tree_data,mut_info) #only merging rows with all data available for now

origin_dpi=0
size=parser$window
slide=parser$window
nec=parser$dpi

Window<-as.data.frame(cbind(Start=seq(origin_dpi,nec-parser$window,slide), End=seq(size,nec,slide)))

mylist.names<-as.character(seq(size,nec,slide))
slid.win<-vector("list", length(seq(origin_dpi,nec-parser$window,slide)))
names(slid.win)<-mylist.names

slid.win = mclapply(1:nrow(Window), function(win) {
  final=mclapply(seq_along(friedrice$origin_dpi), function(branch) {
    if(friedrice$origin_dpi[branch]<=Window$End[win] & 
       friedrice$end_dpi[branch]>=Window$Start[win] ) {
      result = friedrice[branch,]
    } else {
      result=NULL
    }
    return(result)
  }, mc.cores=numCores)
  x=do.call(rbind, compact(final))
  return(data.frame("X"=x))
}, mc.cores=numCores)

#changing names of dataframes to dpi
names(slid.win)<-mylist.names
#changing colnames back to original in dataframes in the list of windows
slid.win = lapply(slid.win, setNames, colnames(friedrice))

#colnames(slid.win$`120`) = colnames(friedrice)  

require(data.table)

#this will put the sequence only in the last window they appear
slid.win.max = rbindlist(slid.win, idcol="window") %>%
  dplyr::mutate(window=as.numeric(window)) %>%
  dplyr::group_by(label) %>%
  dplyr::mutate(max_win=max(window)) %>%
  dplyr::filter(window==max_win)

#this will put the sequence only in the first window they appear
slid.win.min = rbindlist(slid.win, idcol="window") %>%
  dplyr::mutate(window=as.numeric(window)) %>%
  dplyr::group_by(label) %>%
  dplyr::mutate(min_win=min(window)) %>%
  dplyr::filter(window==min_win)

# Gaps often are not isolated; they are found in stretches, so want to consolidate our gapped sites into discrete stretches by
# first identifying mutations that result in gaps (still split up by branch):
deletions.max = filter(slid.win.max, grepl("-$", muts)) %>%
  group_by(label) %>%
  group_split()

deletions.min = filter(slid.win.min, grepl("-$", muts)) %>%
  group_by(label) %>%
  group_split()

for (i in 1) {
  if (is_empty(deletions.max) == TRUE){
    df.max = slid.win.max
    next
  } else {
    deletions.max = lapply(deletions.max, function(x) {
      site_batch = lapply(2:nrow(x), function(y) {
        if(isTRUE(x$site[y] == x$site[y-1]+1)) {
          return(NA)
          } else {
            return(x$site[y])
            }
        })
      site_batch = c(x$site[1], do.call(rbind, site_batch))
      result = cbind(x, site_batch)
      return(result)
      })
# Now categorize gap-related mutations as insertions or deletions
    deletions.max = lapply(deletions.max, function(x) {
      fill(x, site_batch) %>%
        group_by(label, site_batch) %>%
        mutate(batch_length = length(site),
                  batch_end = site_batch + batch_length-1) %>%
      distinct()
    })
# Now take stretches of gaps and concatenate them into one discrete deleted region 
# This is done so that when we look for co-evolving pairs (normally done for individual nucleotide sites),
# We can treat a deletion region as an individual "site" 
    deletions.max = do.call(rbind, lapply(deletions.max, function(x) {
      x <- x %>%
        ungroup() %>%
        mutate(binned = if_else(batch_length>1, 
                              as.character(paste0(site_batch, "-", batch_end)), 
                              as.character(site_batch))) #%>%
    })) %>%
      group_by(binned) %>%
      mutate(num_branches = length(label))
#Now remove indels from the larger mutation datatable and merge it with the deletions table
    require(dplyr)
    df.max = merge(slid.win.max, deletions.max, by = c('label','muts','site','window','max_win','node','parent',
                            'origin_dpi','end_dpi'), all = T)#,'origin_tissue','end_tissue'
write.csv(df.max, paste0(parser$output,"branch-site_mutations_max.win.csv"), quote=F, row.names=F)
  }
}



for (j in 1) {
  if (is_empty(deletions.min) == TRUE){
    df.min = slid.win.min
    next
  } else {
    deletions.min = lapply(deletions.min, function(x) {
      site_batch = lapply(2:nrow(x), function(y) {
        if(isTRUE(x$site[y] == x$site[y-1]+1)) {
          return(NA)
        } else {
          return(x$site[y])
        }
      })
      site_batch = c(x$site[1], do.call(rbind, site_batch))
      result = cbind(x, site_batch)
      return(result)
    })
    # Now categorize gap-related mutations as insertions or deletions
    deletions.min = lapply(deletions.min, function(x) {
      fill(x, site_batch) %>%
        group_by(label, site_batch) %>%
        mutate(batch_length = length(site),
               batch_end = site_batch + batch_length-1) %>%
        distinct()
    })
    # Now take stretches of gaps and concatenate them into one discrete deleted region 
    # This is done so that when we look for co-evolving pairs (normally done for individual nucleotide sites),
    # We can treat a deletion region as an individual "site" 
    deletions.min = do.call(rbind, lapply(deletions.min, function(x) {
      x <- x %>%
        ungroup() %>%
        mutate(binned = if_else(batch_length>1, 
                                as.character(paste0(site_batch, "-", batch_end)), 
                                as.character(site_batch)))
    })) %>%
      group_by(binned) %>%
      mutate(num_branches = length(label))
    #Now remove indels from the larger mutation datatable and merge it with the deletions table
    require(dplyr)
    df.min = merge(slid.win.min, deletions.min, by = c('label','muts','site','window','min_win','node','parent',
                                                       'origin_dpi','end_dpi'), all = T)#,'origin_tissue','end_tissue'
    
    write.csv(df.min, paste0(parser$output,"branch-site_mutations_min.win.csv"), quote=F, row.names=F)
  }
}


# Now, for each branch, lump changing sites into site pairs by using the following function to create all possible pairwise combinations
expand.grid.unique <- function(x, y, include.equals=FALSE) {
  x <- unique(x)
  y <- unique(y)
  g <- function(i)
  {
    z <- setdiff(y, x[seq_len(i-include.equals)])
    if(length(z)) cbind(x[i], z, deparse.level=0)
  }
  do.call(rbind, lapply(seq_along(x), g))
}

# Branches are labelled according to support, but several branches can have the same support value - need to
# Provide an additional layer of naming to distinguish individual branches in downstream analysis:
require(dplyr)
require(bnlearn)

max.per_window = list()
min.per_window = list()

for (g in 1:length(mylist.names)) {
  max.per_window[[g]] = subset(df.max, window == mylist.names[g])
}

for (g in 1:length(mylist.names)) {
  min.per_window[[g]] = subset(df.min, window == mylist.names[g])
}

max.namescol = colnames(max.per_window[[1]])
max.per_window = Reduce(function(...) merge(..., by=max.namescol, all=TRUE), max.per_window)

max.ctable = table(dplyr::select(max.per_window, window, site))

max.ptable = prop.table(max.ctable)

df.max.ctable = as.data.frame.matrix(max.ptable)

max.res = iamb(df.max.ctable)

max.result_pairs = max.res$arcs


min.namescol = colnames(min.per_window[[1]])
min.per_window = Reduce(function(...) merge(..., by=min.namescol, all=TRUE), min.per_window)

min.ctable = table(dplyr::select(min.per_window, window, site))

min.ptable = prop.table(min.ctable)

df.min.ctable = as.data.frame.matrix(min.ptable)

min.res = iamb(df.min.ctable)

min.result_pairs = min.res$arcs

#We can have co-evolving sites where sometimes one site happens without the other
#Grouping into blacklist according to pairs and find pairs with majority having score of zero
#This way, it removes the pairs in which the times the score of 0 is more than the same pair that it was scored as 1
#Gotta remove those pairs that the majority did not have a score of 0 in case there are any

#Used inside the next 2 for loops
removing = function(x) {x - 1}

library(ggplot2)
library(igraph)
library(RColorBrewer)


for (pair in 1) {
  if(is_empty(max.result_pairs)==TRUE){
    print("No pairs of sites initially found in maximum window")
    next
  }else if (nrow(max.ctable)<=1){
    print("All sites have been assigned to the same window. No relevant pairs of sites can be determined")
    next
  } else {
    max.neg_corr=do.call(rbind,mclapply(1:nrow(max.ctable), function(win) {
      do.call(rbind,lapply(1:nrow(max.result_pairs),function(site_pair) {
        val1=max.ctable[win, max.result_pairs[site_pair,1]]
        val2=max.ctable[win, max.result_pairs[site_pair,2]]
        if(val1==0 & val2==0) {
          max.neg_corr = data.frame(from=max.result_pairs[site_pair,1],
                                to=max.result_pairs[site_pair,2], score=0, row.names = NULL)
        } else if (val1 >0 & val2 >0) {
          max.neg_corr = data.frame(from=max.result_pairs[site_pair,1],
                                to=max.result_pairs[site_pair,2], score=1, row.names=NULL)
        } else {
          max.neg_corr = data.frame(from=max.result_pairs[site_pair,1],
                                to=max.result_pairs[site_pair,2], score=0, row.names=NULL)
        }
      }))
    }, mc.cores=numCores))
    
    max.zeros = ddply(max.neg_corr,.(from,to,score),nrow)
    
    indexes.max = as.numeric()
    for (index in 1:(nrow(max.zeros)-1)) {
      index1 = index + 1
      if (max.zeros[index,1]==max.zeros[index1,1] & max.zeros[index,2]==max.zeros[index1,2] & max.zeros[index,4] > max.zeros[index1,4]){
        indexes.max = rbind(indexes.max,index1)
      }
    }
    
    max.zeros= max.zeros %>% slice(-c(indexes.max))
    
    max.minority_zeros = duplicated(max.zeros[,1:2])
    max.minority_zeros_index = which(max.minority_zeros == TRUE)
    
    toremove.max = removing(max.minority_zeros_index)
    max.minority_zeros_index = rbind(max.minority_zeros_index,toremove.max)
    max.zeros = max.zeros %>% slice(-c(max.minority_zeros_index))
    
    #Pairs with majority of score 0
    max.blacklist= max.zeros %>% select(-score,-V1)
    
    # Remove these blacklisted pairs from the original probability table as well
    # After blacklisting, filter pairs based on bootstrap strength:
    max.res=iamb(df.max.ctable, blacklist=max.blacklist)
    bs.max=boot.strength(df.max.ctable, algorithm="iamb")
    
    final.max=do.call(rbind, compact(lapply(split(bs.max, seq(nrow(bs.max))), function(x) {
      x2=paste(x$from,x$to, sep=",")
      bl=do.call(paste, c(max.blacklist, sep=","))
      if(x2 %in% bl) {
        return(NULL)
      } else {
        return(x)
      }
    })) )
    
    # Plots distribution of strengths
    strength_distribution.max = ggplot(final.max, aes(x=strength)) +
      geom_density() + xlab("Strength (posterior)") + ylab("% total pairs of sites") +
      theme_classic()
    ggsave(paste0(parser$output,"_max.win_strength_distribution.png"), strength_distribution.max, width = 8, height = 8, units = "in", dpi = 600)
    
    network.max = strength.plot(max.res, strength = bs.max, main = "Max_window",
                                threshold = 0.50,
                                shape = "circle",
                                layout = "fdp")
    
    network.max1 = graph_from_graphnel(network.max, name = T) 
    network.max1 = delete.edges(network.max1, which(E(network.max1)$weight <0.5))
    png(paste0(parser$output,"_max.win_sites.png"), 1024,1024)
    plot(delete.vertices(simplify(network.max1), degree(network.max1)==0), edge.color="palevioletred4", vertex.label.color="black", vertex.label.font=1,
         vertex.color = "lightpink2", edge.arrow.size = 1, vertex.size = 30,
         vertex.label.cex=2, edge.lty = 1)
    dev.off()
    
  }
}

for (pair in 1) {
  if(is_empty(min.result_pairs)==TRUE){
    print("No pairs of sites found in minimum window")
    next
  }else if (nrow(min.ctable)<=1){
    print("No relevant pairs of sites can be determined")
    next
  } else {  
    min.neg_corr=do.call(rbind,mclapply(1:nrow(min.ctable), function(win) {
      do.call(rbind,lapply(1:nrow(min.result_pairs),function(site_pair) {
        val1=min.ctable[win, min.result_pairs[site_pair,1]]
        val2=min.ctable[win, min.result_pairs[site_pair,2]]
        if(val1==0 & val2==0) {
          min.neg_corr = data.frame(from=min.result_pairs[site_pair,1],
                                to=min.result_pairs[site_pair,2], score=0, row.names = NULL)
        } else if (val1 >0 & val2 >0) {
          min.neg_corr = data.frame(from=min.result_pairs[site_pair,1],
                                to=min.result_pairs[site_pair,2], score=1, row.names=NULL)
        } else {
          min.neg_corr = data.frame(from=min.result_pairs[site_pair,1],
                                to=min.result_pairs[site_pair,2], score=0, row.names=NULL)
        }
      }))
    }, mc.cores=numCores))
    
    min.zeros = ddply(min.neg_corr,.(from,to,score),nrow)
    
    indexes.min = as.numeric()
    for (index in 1:(nrow(min.zeros)-1)) {
      index1 = index + 1
      if (min.zeros[index,1]==min.zeros[index1,1] & min.zeros[index,2]==min.zeros[index1,2] & min.zeros[index,4] > min.zeros[index1,4]){
        indexes.min = rbind(indexes.min,index1)
      }
    }
    
    min.zeros= min.zeros %>% slice(-c(indexes.min))
    
    min.minority_zeros = duplicated(min.zeros[,1:2])
    min.minority_zeros_index = which(min.minority_zeros == TRUE)
    
    toremove.min = removing(min.minority_zeros_index)
    min.minority_zeros_index = rbind(min.minority_zeros_index,toremove.min)
    min.zeros = min.zeros %>% slice(-c(min.minority_zeros_index))
    
    #Pairs with majority of score 0
    min.blacklist=min.zeros %>% select(-score,-V1)
    
    # Remove these blacklisted pairs from the original probability table as well
    # After blacklisting, filter pairs based on bootstrap strength:
    min.res=iamb(df.min.ctable, blacklist=min.blacklist)
    bs.min=boot.strength(df.min.ctable, algorithm="iamb")
    
    final.min=do.call(rbind, compact(lapply(split(bs.min, seq(nrow(bs.min))), function(x) {
      x2=paste(x$from,x$to, sep=",")
      bl=do.call(paste, c(min.blacklist, sep=","))
      if(x2 %in% bl) {
        return(NULL)
      } else {
        return(x)
      }
    })) )
    
    # Plots distribution of strengths
    strength_distribution.min = ggplot(final.min, aes(x=strength)) +
      geom_density() + xlab("Strength (posterior)") + ylab("% total pairs of sites") +
      theme_classic()
    ggsave(paste0(parser$output,"_min.win_strength_distribution.png"), strength_distribution.min, width = 8, height = 8, units = "in", dpi = 600)
   
    network.min = strength.plot(min.res, strength = bs.min, main = "Min_window",
                                threshold = 0.50,
                                shape = "circle",
                                layout = "fdp") 
    network.min1 = graph_from_graphnel(network.min, name = T) 
    network.min1 = delete.edges(network.min1, which(E(network.min1)$weight <0.5))
    png(paste0(parser$output,"_min.win_sites.png"), 1024,1024)
    plot(delete.vertices(simplify(network.min1), degree(network.min1)==0), edge.color="darkblue", vertex.label.color="black", vertex.label.font=1,
         vertex.color = "lightblue", edge.arrow.size = 1, vertex.size = 30,
         vertex.label.cex=2, edge.lty = 1)
    dev.off()
    
  }
}

