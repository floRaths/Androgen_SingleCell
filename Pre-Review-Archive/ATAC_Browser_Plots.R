function(){
options(stringsAsFactors=FALSE)

def.args = c("/avicenna/mkarimzadeh/cisTransBreastProject/splittedBams/footPrintAnalysis/20201013/plots",
             "/avicenna/mkarimzadeh/cisTransBreastProject/archrRun/GenesAndTheirEnhancerPositions_allCorrelations.RDS",
             "/avicenna/mkarimzadeh/cisTransBreastProject/splittedBams/footPrintAnalysis/20201013/plots/geneAnnotations.RDS",
             "/avicenna/mkarimzadeh/cisTransBreastProject/splittedBams/footPrintAnalysis/20201013/plots/exonAnnotations.RDS",
             "/avicenna/mkarimzadeh/cisTransBreastProject/splittedBams/peaks/LUM_HR-pos_Trans/narrowPeaks/LUM_HR-pos_Trans_treat_pileup.bigWig",
             "/avicenna/mkarimzadeh/cisTransBreastProject/splittedBams/peaks/LUM_HR-pos_Cis/narrowPeaks/LUM_HR-pos_Cis_treat_pileup.bigWig",
             "/avicenna/mkarimzadeh/cisTransBreastProject/splittedBams/footPrintAnalysis/20201013/plots/motifPositions.cisBP.RDS",
             "--genes", "AR", "PGR",
             "--motif", "AR", "JUND")

args = commandArgs(TrailingOnly=TRUE)


if(length(args) < 2){
  cat("Usage: Rscript plot_peaks_and_motifs.R <outDir> <peakInfo> <geneTable> <exonTable> <bwtrans> <bwcis> <motifPoses> --genes <...> <...> --motifs <...> <...>\n")
  cat("Using default arguments\n")
  print(def.args)
  args = def.args
}

}

col <- pal_d3("category20")(20)

# plotting output function
create_plot = function(P, path, width, height){
  #png(sprintf("%s.png", path), res=200, width=width*200, height=height*200)
  #plot(P)
  #dev.off()
  pdf(sprintf("%s.pdf", path), width=width, height=height)
  plot(P)
  dev.off()
}
# moving average function
ma <- function(arr, n = 15){
  res = arr
  for(i in n:length(arr)){
    res[i] = mean(arr[(i-n):i])
  }
  res
}
make_smooth_peak_signal = function(bwobj1, bwobj2, p2gdf, window = 10){
  tempobj = GRanges(seqnames=p2gdf$ChromAtac,
                    ranges=IRanges(start=p2gdf$StartAtac - 1000,
                                   end=p2gdf$EndAtac + 1000))
  bw1 = subsetByOverlaps(bwobj1, tempobj)
  bw2 = subsetByOverlaps(bwobj2, tempobj)
  print("Subsetted initial input")
  list.bwobjs = list(Cis=bw1, Trans=bw2)
  idxs_atac = unique(p2gdf[,"idxATAC"])
  list.data = lapply(idxs_atac, function(idx){
    print(idx)
    start = p2gdf$StartAtac[p2gdf$idxATAC==idx][1] - 1000
    end = p2gdf$EndAtac[p2gdf$idxATAC==idx][1] + 1000
    tempobj = GRanges(seqnames=p2gdf$ChromAtac[p2gdf$idxATAC==idx][1],
                      ranges=IRanges(
                        start=seq(start, end),
                        end=seq(start, end)+1))
    list.ad = list()
    for(i in 1:2){
      adname = names(list.bwobjs)[i]
      bwobj = list.bwobjs[[i]]
      idxoverlap = as.data.frame(findOverlaps(bwobj, tempobj))
      tempdf = as.data.frame(tempobj)
      tempdf$score = 0
      tempdf$score[idxoverlap[,2]] = bwobj$score[idxoverlap[,1]]
      tempdf$idxATAC = idx
      tempdf$Group = adname
      tempdf$smoothedSignal = ma(tempdf$score, window * 5)
      tempdf$Resolution = "1 bp"
      # Now use every 10 bp
      tempobj2 = GRanges(seqnames=p2gdf$ChromAtac[p2gdf$idxATAC==idx][1],
                         ranges=IRanges(start=seq(start, end, by=window),
                                        end=seq(start, end, by=window)+window + 1))
      idxoverlap = as.data.frame(findOverlaps(tempobj2, tempobj))
      newdf = as.data.frame(tempobj2)
      for(variable in c("score", "smoothedSignal")){
        newdf[,variable] = sapply(unique(idxoverlap[,1]), function(i){
          return(mean(tempdf[idxoverlap[idxoverlap[,1]==i,2],variable]))})
      }
      for(variable in c("idxATAC", "Group")){
        newdf[,variable] = unique(tempdf[,variable])
      }
      newdf$Resolution = paste(window, "bp")
      newdf2 = rbind(tempdf, newdf)
      list.ad = c(list.ad, list(newdf2))
    }
    addf = do.call("rbind", list.ad)
    return(addf)
  })
  outdf = do.call("rbind", list.data)
  return(outdf)
}

#make_smooth_peak_signal = function(bwobj1, bwobj2, p2gdf, window=10, START=NA, END=NA){
  if(is.na(START)){
    start = min(p2gdf$StartAtac) - 1000
    end = max(p2gdf$EndAtac) + 1000
  }else{
    start = START
    end = END
  }
  tempobj = GRanges(seqnames=p2gdf$ChromAtac[1],
                    ranges=IRanges(start=seq(start, end),
                                   end=seq(start, end) + 1))
  bw1 = subsetByOverlaps(bwobj1, tempobj)
  bw2 = subsetByOverlaps(bwobj2, tempobj)
  print("Subsetted initial input")
  list.bwobjs = list(Cis=bw1, Trans=bw2)
  idx = NA
  list.ad = list()
  for(i in 1:2){
    adname = names(list.bwobjs)[i]
    bwobj = list.bwobjs[[i]]
    idxoverlap = as.data.frame(findOverlaps(bwobj, tempobj))
    tempdf = as.data.frame(tempobj)
    tempdf$score = 0
    tempdf$score[idxoverlap[,2]] = bwobj$score[idxoverlap[,1]]
    tempdf$idxATAC = NA
    tempdf$Group = adname
    tempdf$smoothedSignal = ma(tempdf$score, window * 5)
    tempdf$Resolution = "1 bp"
    # Now use every 10 bp
    tempobj2 = GRanges(seqnames=p2gdf$ChromAtac[1],
                       ranges=IRanges(start=seq(start, end, by=window),
                                      end=seq(start, end, by=window)+window + 1))
    
    idxoverlap = as.data.frame(findOverlaps(tempobj2, tempobj))
    newdf = as.data.frame(tempobj2)
    for(variable in c("score", "smoothedSignal")){
      newdf[,variable] = sapply(unique(idxoverlap[,1]), function(i){
        return(mean(tempdf[idxoverlap[idxoverlap[,1]==i,2],variable]))})
    }
    for(variable in c("idxATAC", "Group")){
      newdf[,variable] = unique(tempdf[,variable])
    }
    newdf$Resolution = paste(window, "bp")
    # newdf2 = rbind(tempdf, newdf)
    list.ad = c(list.ad, list(newdf))
  }
  addf = do.call("rbind", list.ad)
  idxs_atac = unique(p2gdf[,"idxATAC"])
  list.data = lapply(idxs_atac, function(idx){
    print(idx)
    start = p2gdf$StartAtac[p2gdf$idxATAC==idx][1] - 1000
    end = p2gdf$EndAtac[p2gdf$idxATAC==idx][1] + 1000
    tempobj = GRanges(seqnames=p2gdf$ChromAtac[p2gdf$idxATAC==idx][1],
                      ranges=IRanges(
                        start=seq(start, end),
                        end=seq(start, end)+1))
    list.ad = list()
    for(i in 1:2){
      adname = names(list.bwobjs)[i]
      bwobj = list.bwobjs[[i]]
      idxoverlap = as.data.frame(findOverlaps(bwobj, tempobj))
      tempdf = as.data.frame(tempobj)
      tempdf$score = 0
      tempdf$score[idxoverlap[,2]] = bwobj$score[idxoverlap[,1]]
      tempdf$idxATAC = idx
      tempdf$Group = adname
      tempdf$smoothedSignal = ma(tempdf$score, window * 5)
      tempdf$Resolution = "1 bp"
      # Now use every 10 bp
      tempobj2 = GRanges(seqnames=p2gdf$ChromAtac[p2gdf$idxATAC==idx][1],
                         ranges=IRanges(start=seq(start, end, by=window),
                                        end=seq(start, end, by=window)+window + 1))
      idxoverlap = as.data.frame(findOverlaps(tempobj2, tempobj))
      newdf = as.data.frame(tempobj2)
      for(variable in c("score", "smoothedSignal")){
        newdf[,variable] = sapply(unique(idxoverlap[,1]), function(i){
          return(mean(tempdf[idxoverlap[idxoverlap[,1]==i,2],variable]))})
      }
      for(variable in c("idxATAC", "Group")){
        newdf[,variable] = unique(tempdf[,variable])
      }
      newdf$Resolution = paste(window, "bp")
      newdf2 = rbind(tempdf, newdf)
      list.ad = c(list.ad, list(newdf2))
    }
    addf = do.call("rbind", list.ad)
    return(addf)
  })
  outdf = do.call("rbind", list.data)
  outdf = rbind(addf, outdf)
  return(outdf)
}

require(parallel)
require(GRanges)
require(rtracklayer)
require(ggplot2)
require(reshape2)
require(scales)
require(cowplot)
require(ggrepel)
require(RColorBrewer)

celltype = "Adipocyte"

args = c("~/Documents/ATAC_Data/Output/plots/",
         "~/Documents/ATAC_Data/scriptData/GenesAndTheirEnhancerPositions_allCorrelations.RDS",
         "~/Documents/ATAC_Data/scriptData/geneAnnotations.RDS",
         "~/Documents/ATAC_Data/scriptData/exonAnnotations.RDS",
         paste0("~/Documents/ATAC_Data/bigWigFiles/", celltype, "_Trans_treat_pileup.bigWig"),
         paste0("~/Documents/ATAC_Data/bigWigFiles/", celltype, "_Cis_treat_pileup.bigWig"),
         "~/Documents/ATAC_Data/scriptData/motifPositions.cisBP.RDS",
         "--genes", "FMO5",
         "--motif", "AR", "PGR", "ESR1", "NR4A1", "EGR1", "NR3C1")


outdir         = args[1]
p2gdf          = readRDS(args[2])
geneobj        = readRDS(args[3])
exonobj        = readRDS(args[4])

cat(sprintf("Loading %s\n", args[5]))
bwtrans        = import(args[5])
cat(sprintf("Loading %s\n", args[6]))
bwcis          = import(args[6])

motifPositions = readRDS(args[7])
idx_motif      = which(args == "--motif")
#genes          = args[9:(idx_motif - 1)]
genes          = c("MLXIPL", "NR3C1")
motifs         = args[(idx_motif + 1):length(args)]


## motif selection
MOTIFS = unique(unlist(lapply(motifs, function(motif){
  admotifs = names(motifPositions)[grep(motif, names(motifPositions))]
  return(admotifs)})))

cat(sprintf("Will use the following motifs: \n %s\n", paste(MOTIFS, collapse="\n")))
motifs = MOTIFS


refseqdf = as.data.frame(geneobj)
exondf   = as.data.frame(exonobj)
p2gdf_celltype = p2gdf[p2gdf$GeneRna %in% genes, ]

### takes_time
sigdf        = make_smooth_peak_signal(bwcis, bwtrans, p2gdf_celltype)#, START = min(p2gdf_celltype$StartAtac - 25000), END = max(p2gdf_celltype$StartAtac) + 25000)

idx_singlebp = which(sigdf$Resolution == "1 bp")
sigdf_1bp    = sigdf[idx_singlebp, ]
sigdf_10bp   = sigdf[idx_singlebp, ]

# Granges object for easier access
sigobj       = GRanges(seqnames=sigdf_1bp$seqnames,
                 ranges=IRanges(start=sigdf_1bp$start,
                                end=sigdf_1bp$end))


list.motif.sig = mclapply(motifs, function(motif){
  print(sprintf("Adding %s", motif))
  motifobj = motifPositions[[motif]]
  tempdf_motif = as.data.frame(motifobj)
  motifobj_temp = GRanges(seqname=tempdf_motif$seqnames,
                          ranges=IRanges(start=motifobj@ranges@start - 150,
                                         end=motifobj@ranges@start + motifobj@ranges@width + 150))
  idxdf = as.data.frame(findOverlaps(sigobj, motifobj_temp))
  if(nrow(idxdf) > 0){
    sigdf_motif = sigdf_1bp[idxdf[,1], ]
    sigdf_motif$MotifName = motif
    sigdf_motif$MotifCenter = tempdf_motif[idxdf[,2],2] + round(tempdf_motif[idxdf[,2],4] / 2)
    sigdf_motif$RelativeStart = sigdf_motif$start - sigdf_motif$MotifCenter
    sigdf_motif$MotifSignal = tempdf_motif[idxdf[,2],"score"]
    return(sigdf_motif)
  }else{
    return(NULL)
  }
}, mc.cores = 16)

motifsigdf = do.call("rbind", list.motif.sig)


motifsigdf$GeneRna = sapply(motifsigdf$idxATAC, function(x){
  return(p2gdf_celltype$GeneRna[p2gdf_celltype$idxATAC==x][1])})



#for(gene in unique(p2gdf_celltype$GeneRna)){
funtion(gene, display_mots) {  
  
  idxatacs = p2gdf_celltype[p2gdf_celltype$GeneRna==gene, "idxATAC"]
  
  p2gdf_gene = p2gdf_celltype[p2gdf_celltype$GeneRna==gene,]
  p2gdf_test = p2gdf_gene[1,]
  p2gdf_test$StartAtac = min(p2gdf_gene$StartAtac - 5000)
  p2gdf_test$EndAtac = max(p2gdf_gene$EndAtac + 5000)
  
  #takes the most time
  sigdf_cur_gene = make_smooth_peak_signal(bwcis, bwtrans, p2gdf_test, window = 50)#, 
                                           #min(p2gdf$StartAtac - 25000), 
                                           #max(p2gdf$StartAtac) + 25000)
  
  
  # Process gene df
  sigdf_gene = sigdf_10bp[sigdf_10bp$idxATAC %in% idxatacs, ]
  sigdf_gene = sigdf_gene[sigdf_gene$smoothedSignal!=0, ]
  sigdf_gene$minSignal = ifelse(grepl("Cis", sigdf_gene$Group), -1 * sigdf_gene$smoothedSignal, 0)
  sigdf_gene$maxSignal = ifelse(grepl("Cis", sigdf_gene$Group), 0, sigdf_gene$smoothedSignal)
  sigdf_gene$Label = "Accessibility"
  maxval = max(abs(sigdf_gene$smoothedSignal))
  
  # Process motifs
  sigdf_gene_motif = motifsigdf[motifsigdf$idxATAC %in% idxatacs, ]
  motifscores = sapply(unique(sigdf_gene_motif$MotifName), function(x){
    return(max(sigdf_gene_motif$smoothedSignal[sigdf_gene_motif$MotifName==x]))})
  
  # motifscores = motifscores[grep("(ZBTB|BPTF|HLF|CEBP)", names(motifscores))]
  if(length(motifscores) <= 9){
    use_motifs = names(motifscores)
  }else{
    use_motifs = names(motifscores)[order(motifscores, decreasing=TRUE)][1:8]
  }
  
  #use_motifs = union(use_motifs, "BPTF_812")
  #use_motifs = c("CUX2_292", "AR_689", "NR3C1_666", "FOXA1_357")
  #use_motifs = c("NFIC_740", "NR3C1_666", "TP63_704")
  #use_motifs = unique %>% unique()
  use_motifs = c("AR_689", "NR4A1_671", "BPTF_812",   "EGR1_195")
  
  sigdf_gene_motif_1 = sigdf_gene_motif[sigdf_gene_motif$MotifName %in% use_motifs, ]
  sigdf_gene_motif_1$MotifName = factor(sigdf_gene_motif_1$MotifName, levels=use_motifs)
  sigdf_gene_motif_1$Grouping = paste(sigdf_gene_motif_1$Group, sigdf_gene_motif_1$MotifCenter, sigdf_gene_motif_1$idxATAC)
  textdf = unique(sigdf_gene_motif_1[,c("Group", "MotifCenter", "idxATAC", "MotifName", "Grouping")])
  textdf$RelativeStart = seq(-100, to=100, length.out=nrow(textdf))
  
  textdf$smoothedSignal = sapply(textdf$Grouping, function(x){
    cur_st = textdf$RelativeStart[textdf$Grouping==x]
    tempdf = sigdf_gene_motif_1[sigdf_gene_motif_1$Grouping==x, ]
    sim_val = tempdf[which.min(abs(tempdf$RelativeStart - cur_st)), "smoothedSignal"][1]
    return(sim_val)})
  
  textdf$MotifName = factor(textdf$MotifName, levels=use_motifs)
  gene_start = refseqdf[!is.na(refseqdf$symbol) & refseqdf$symbol==gene ,2][1]
  gene_end = refseqdf[!is.na(refseqdf$symbol) & refseqdf$symbol==gene ,3][1]
  exons_gene = exondf[!is.na(exondf$symbol) & exondf$symbol==gene, ]
  cur_genedf = refseqdf[!is.na(refseqdf$symbol) & refseqdf$symbol==gene , ]
  chevrondf = data.frame(start=seq(gene_start, gene_end, length.out=50),
                         end=seq(gene_start, gene_end, length.out=50) + 1,
                         score=-maxval * 0.75,
                         Character=ifelse(cur_genedf$strand=="+", ">", "<"))

  
  sigdf_cur_gene = sigdf_cur_gene[sigdf_cur_gene$Resolution != "1 bp", ]
  sigdf_cur_gene$minSignal = ifelse(grepl("Cis", sigdf_cur_gene$Group), -1 * sigdf_cur_gene$smoothedSignal, 0)
  sigdf_cur_gene$maxSignal = ifelse(grepl("Cis", sigdf_cur_gene$Group), 0, sigdf_cur_gene$smoothedSignal)
  
#}
  


  #prefix = 189
  mod_maxval = (maxval - 0)
  
  
  # Mark gene promoter
  #p = 
    ggplot() +
    #geom_rect(data=sigdf_gene, aes(xmin=start, xmax=end), ymin=-mod_maxval, ymax=(-mod_maxval * 0.95), fill="grey40") +
    geom_rect (data=cur_genedf, aes(xmin=start, xmax=end), ymin=(-mod_maxval * 0.8), ymax=(-mod_maxval * 0.7),   fill="deepskyblue4", alpha = 0.5) + #gene body
    #geom_rect (data=exons_gene, aes(xmin=start, xmax=end), ymin=(-mod_maxval * 0.85), ymax=(-mod_maxval * 0.65), fill="deepskyblue4", alpha = 0.25) + #exons
    geom_rect (data=sigdf_cur_gene, aes(xmin=start, xmax=end, ymin=minSignal, ymax=maxSignal, fill=Group)) +
    geom_segment(data=unique(top_n(group_by(select(textdf, idxATAC, MotifCenter, MotifName), idxATAC, MotifName), 1, MotifCenter)), 
                 aes(x = MotifCenter, y = mod_maxval, xend = MotifCenter, yend = (0.85*mod_maxval), colour = MotifName)) +
    
    annotate("text", x=chevrondf$start, y=(-mod_maxval * 0.75), colour=rep("grey30", nrow(chevrondf)),
             label=chevrondf$Character) +
    
    scale_fill_manual(name="", values=c("#86007D", "#FFA52C")) +
    scale_color_manual(values = col) +
    theme_bw(base_size=18) +
    
    scale_y_continuous(name=paste(gene, as.character(cur_genedf$seqnames)), limits=c(-mod_maxval, mod_maxval), expand = c(0,0)) +
    #scale_x_continuous(name=paste(gene, "Genomic position", as.character(cur_genedf$seqnames)), labels=comma) + 
                       #limits=c(as.numeric(paste0(prefix, 125, "000")), as.numeric(paste0(prefix, 260, "000")))
                       #) +
      
    #coord_cartesian(xlim = c(as.numeric(paste0("112", "800", "000")), 
    #                         as.numeric(paste0("113", "160" , "000")))) +
      
    theme(text = element_text(family = "Lato"),
          axis.title.x = element_blank(), 
          legend.position = "left",
          #legend.position = "none"
          )
  
  
}
    
#p <- FMO5 + SORD + ANPEP + ABCC11 + IQGAP2 + PPIF + patchwork::plot_layout(ncol = 1)    

p %>% save_x(data = ., name = "TCF7_AR_Adipocytes", 1.5, 10, 3, svg = T) 
  

#### footprints by peak  
  
  for(each_peak in intersect(sigdf_gene_motif$idxATAC, sigdf_gene$idxATAC)){
    
    sigdf_gene_motif_1 = sigdf_gene_motif[sigdf_gene_motif$idxATAC == each_peak, ]
    ad_pos = paste(sigdf[1,1], min(sigdf[1,2]), max(sigdf[1,3]), sep=".")
    outpath = sprintf("%s/%s_peakFootPrint_peak_%d_at_%s",
                      outdir, gene,  each_peak, ad_pos)
    sigdf_gene_1 = sigdf_gene[sigdf_gene$idxATAC==each_peak, ]
    motifscores = sapply(unique(sigdf_gene_motif_1$MotifName), function(x){
      return(max(sigdf_gene_motif$smoothedSignal[sigdf_gene_motif_1$MotifName==x]))})
    if(length(motifscores) <= 8){
      use_motifs = names(motifscores)
    }else{
      use_motifs = names(motifscores)[order(motifscores, decreasing=TRUE)][1:8]
    }
    use_motifs = union(use_motifs)
    #use_motifs = c("CUX2_292", "AR_689", "NR3C1_666", "FOXA1_357", "JUNB_139")
    use_motifs = filter
    
    sigdf_gene_motif_1 = sigdf_gene_motif_1[sigdf_gene_motif_1$MotifName %in% use_motifs, ]
    sigdf_gene_motif_1$MotifName = factor(sigdf_gene_motif_1$MotifName, levels=use_motifs)
    sigdf_gene_motif_1$Grouping = paste(sigdf_gene_motif_1$Group, sigdf_gene_motif_1$MotifCenter, sigdf_gene_motif_1$idxATAC)
    textdf = unique(sigdf_gene_motif_1[,c("Group", "MotifCenter", "idxATAC", "MotifName", "Grouping")])
    textdf$RelativeStart = seq(-100, to=100, length.out=nrow(textdf))
    textdf$smoothedSignal = sapply(textdf$Grouping, function(x){
      cur_st = textdf$RelativeStart[textdf$Grouping==x]
      tempdf = sigdf_gene_motif_1[sigdf_gene_motif_1$Grouping==x, ]
      sim_val = tempdf[which.min(abs(tempdf$RelativeStart - cur_st)), "smoothedSignal"][1]
      return(sim_val)})
    textdf$MotifName = factor(textdf$MotifName, levels=use_motifs)
    ## Add annotation of the motifs to sigdf_gene
    
    P2 = 
      ggplot(sigdf_gene_motif_1, aes(x=RelativeStart, y=smoothedSignal, colour=MotifName, group=Grouping)) +
      geom_line(alpha=0.5) +
      scale_colour_brewer(palette="Set1") +
      # geom_point(alpha=0.5) +
      geom_text_repel(data=textdf, aes(label=MotifName), force=10) +
      theme_bw(base_size=18) +
      ggtitle(sprintf("Peaks of %s", gene)) +
      facet_grid(~Group, scales="free") +
      theme(legend.position="bottom") +
      scale_x_continuous("Relative genomic position to peak center", limits=c(-150, 150)) +
      ylab("Pseudo-bulk ATAC-seq FPM")
    
    P = 
      ggplot(sigdf_gene[sigdf_gene$idxATAC %in% sigdf_gene_motif_1$idxATAC,], aes(xmin=start, xmax=end, ymin=minSignal, ymax=maxSignal, fill=Group)) +
      geom_rect() +
      geom_vline(data=textdf, aes(xintercept=MotifCenter, colour=MotifName)) +
      scale_fill_manual(name="", values=c("#86007D", "#FFA52C")) +
      scale_colour_brewer(palette="Set1") +
      facet_grid(Label~idxATAC) +
      theme_bw(base_size=18) +
      scale_y_continuous(name="Chromatin accessibility", limits=c(-maxval, maxval)) +
      scale_x_continuous(name="Genomic position", labels=comma) +
      theme(legend.position="none",
            axis.text.x=element_text(angle=45, hjust=1))
    
      
      P3 = plot_grid(P2, P, ncol=1, rel_heights=c(1, 2), align="hv")
    
      create_plot(P3, outpath,
                width=10, height=10)
  }
}

each_peak


P2 + P + patchwork::plot_layout(ncol = 2, widths = c(1,2))




p <- 
PGR_2 + theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank()) + 
PGR_3 + theme(legend.position = "none", axis.title = element_blank()) + 
CUX2_2 + theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank()) + 
CUX2_3 + theme(legend.position = "none", axis.title = element_blank()) + 
P2 + theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank()) + 
  P + theme(legend.position = "none", axis.title = element_blank()) + 
  
  patchwork::plot_layout(ncol = 2, widths = c(1,2))
