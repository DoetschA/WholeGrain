#' create barplot of most abundant taxa
#' 
#' This script uses functions of phyloseq and ggplot to display bargraphs showing the relative abundance of taxa in a phyloseq object.
#' 
#' @param physeq A phyloseq object containg at least abundance (otu_table) and taxonomy (tax_table) data. Metadata (sample_data) are required to be able to group samples.
#' @param k Display the k most abundant taxa in the specified dataset (default: 17).
#' @param usePalette Specifies the color palette used to fill the bar sections representing taxa (Default: a 17-color scheme built with RColorBrewer - c(brewer.pal(9,"Set1"), brewer.pal(8,"Dark2")) ).
#' @param taxlevel Taxa (OTUs) are grouped at this taxonomic level for the purpose of coloring. Individual OTUs will still be shown as separate sections of the bar filled with the same color, unless colorfill = T (Default: "Genus").
#' @param relative If T shows relative abundances summing up to 100 % (although only the top k taxa will be shown). If F the bars will be the absolute counts as specified in otu_table (Default: T).
#' @param groupBy Group samples by this factor present as sample_variable of the phyloseq object.
#' @param colorfill If T uses the same color for fill and border of the bars, effectively removes the visibility of single OTUs and only shows the abundance of the taxa at the specified taxlevel (Default: F).
#' @examples 
#' bargraph_top_taxa(physeq, k = 17, taxlevel = "Genus") # show the relative abundance of the top 17 genera of all samples
#' bargraph_top_taxa(physeq, k = 17, taxlevel = "Genus", groupBy = "group") # show the average relative abundance of the top 17 genera of different sample 'group's
#' 
#' 
bargraph_top_taxa <- function(physeq, k = 17, usePalette = c(brewer.pal(9,"Set1"), brewer.pal(8,"Dark2")), taxlevel = "Genus", relative = T, groupBy = NULL, colorfill = F){
  require(phyloseq)
  require(RColorBrewer)
  require(ggplot2)
  
  # transform absolute to relative counts (recommended for grouped samples to avoid  means weighted by group size)
  if(relative){
    physeq <- transform_sample_counts(physeq, function(x) 100 * x/sum(x))
  }
  
  # group by a variable and determine mean abundances
  if(!is.null(groupBy)){
    groupLevels <- levels(sample_data(physeq)@.Data[[which(names(sample_data(physeq)) == groupBy)]])
    # merge samples and replace the summed abundances by summed abundances (mean is not implemented correctly!)
    physeq <- merge_samples(physeq, groupBy)
    # repeating the transformation gets the values back to percent, effectively calculating the mean
    if(relative){
      physeq <- transform_sample_counts(physeq, function(x) 100 * x/sum(x))
    }
  }
  
  
  # find top k unique taxa on requested level, or if k is not numeric use the provided list of taxa names
  list_taxa <- as.character(unique(tax_table(physeq)[names(sort(taxa_sums(physeq), decreasing = T)), taxlevel]))
  if(is.numeric(k)){
    # numeric k => use top k taxa
    if(length(list_taxa) >= k){
      keep_taxa <- list_taxa[1:k]
    } else {
      keep_taxa <- list_taxa
    }
  } else {
    keep_taxa <- list_taxa[list_taxa %in% k]
  }

  
  physeq_selection <- prune_taxa(tax_table(physeq)[,taxlevel] %in% keep_taxa, physeq)
  
  # barplot of top k taxa  
  if(relative) { ytext <- "relative abundance (%)" } else { ytext <- "read counts"}
  
  if(is.null(groupBy)){
    p_bar <- plot_bar(physeq_selection, fill=taxlevel) + ylab(ytext) + theme(text = element_text(size=14)) + scale_fill_manual(values = usePalette) + scale_color_manual(values = usePalette)
  } else {
    p_bar <- plot_bar(physeq_selection, fill=taxlevel) + ylab(ytext) + theme(text = element_text(size=14)) + scale_fill_manual(values = usePalette) + scale_color_manual(values = usePalette) + scale_x_discrete(limits = groupLevels)
  }
  
  if(colorfill){
    eval(parse(text=paste("p_bar + geom_bar(aes(color=", taxlevel, ", fill=", taxlevel, "), stat='identity', position='stack')", sep="")))
  } else {
    p_bar
  }
  
}