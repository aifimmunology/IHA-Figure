# ##### FUNCTIONS ######


setCutoffsMarkers <- function(dt=dat, 
                              pa=images.path) {
  # Identify cutoffs for and label CXCR3+ and CD69+ B cells
  # based on biaxial expression plot evaluation
  # Inputs:
  #   dt - data.table
  #   pa - path to images folder
  # Outputs:
  #   dt - data.table
  
  get_density <- function(x, y, ...) {
    dens <- MASS::kde2d(x, y, ...)
    ix <- findInterval(x, dens$x)
    iy <- findInterval(y, dens$y)
    ii <- cbind(ix, iy)
    return(dens$z[ii])
  }  
  
  print("Setting cutoffs CXCR3+ B cells")
  
  min(table(dt$timepoint))
  n.sample <- 3000
  set.seed(888)
  subsampled <- dt[non_B==F, .SD[sample(.N, n.sample)], by=.(timepoint)] 
  
  ggplot(subsampled, aes(CD183, HLA_DR, fill=CD19)) + 
    geom_point(color="black", pch=21) + 
    scale_fill_viridis(option="B") + 
    xlim(0.4, 1.1) +
    theme(legend.position="right") + 
    theme_minimal() + 
    geom_vline(xintercept=0.85)
  
  ggplot(subsampled, aes(CD183, HLA_DR, fill=timepoint)) + 
    geom_point(color="black", pch=21) + 
    xlim(0.4, 1) +
    theme(legend.position="right") + 
    theme_minimal() + 
    geom_vline(xintercept=0.85)
  
  dt[, `:=`(CXCR3_pos=F)]
  
  dt[CD183>0.85, CXCR3_pos:=T]
  
  dt[, .N, by=.(CXCR3_pos)]
  
  n.sample <- 3000
  set.seed(888)
  subsampled <- dt[non_B==F, .SD[sample(.N, n.sample)], by=.(timepoint)]
  
  ggplot(subsampled, aes(CD183, HLA_DR, fill=timepoint)) + 
    geom_point(color="black", pch=21) + 
    xlim(0.4, 1) +
    theme(legend.position="right") + 
    theme_minimal() + 
    geom_vline(xintercept=0.85)
  ggsave(paste0(pa, "CXCR3_pos_B_cell_biaxial.eps"), width=8, height=6)
  ggsave(paste0(pa, "CXCR3_pos_B_cell_biaxial.png"), width=8, height=6)
  
  print("Setting cutoffs CD69+ B cells")
  
  min(table(dt$timepoint))
  n.sample <- 1000
  set.seed(888)
  subsampled <- dt[non_B==F, .SD[sample(.N, n.sample)], by=.(timepoint)] 
  
  ggplot(subsampled, aes(CD69, CD19, fill=HLA_DR)) + 
    geom_point(color="black", pch=21) + 
    scale_fill_viridis(option="B") + 
    xlim(0, 1) +
    theme(legend.position="right") + 
    theme_minimal() + 
    geom_vline(xintercept=0.45)
  
  ggplot(subsampled, aes(CD69, HLA_DR, fill=timepoint)) + 
    geom_point(color="black", pch=21) + 
    xlim(0, 1) +
    theme(legend.position="right") + 
    theme_minimal() + 
    geom_vline(xintercept=0.45)
  
  dt[, `:=`(CD69_pos=F)]
  
  dt[CD69>0.45, CD69_pos:=T]
  
  dt[, .N, by=.(CD69_pos)]
  
  n.sample <- 1000
  set.seed(888)
  subsampled <- dt[non_B==F, .SD[sample(.N, n.sample)], by=.(timepoint)]
  
  ggplot(subsampled, aes(CD69, HLA_DR, fill=timepoint)) + 
    geom_point(color="black", pch=21) + 
    xlim(0, 1) +
    theme(legend.position="right") + 
    theme_minimal() + 
    geom_vline(xintercept=0.45)
  ggsave(paste0(pa, "CD69_pos_B_cell_biaxial.eps"), width=8, height=6)
  ggsave(paste0(pa, "CD69_pos_B_cell_biaxial.png"), width=8, height=6)
  
  
  print("Setting cutoffs CD27+ B cells")
  
  min(table(dt$timepoint))
  n.sample <- 10000
  set.seed(888)
  subsampled <- dt[non_B==F, .SD[sample(.N, n.sample)], by=.(timepoint)] 
  
  ggplot(subsampled, aes(CD27, CD19, fill=HLA_DR)) + 
    geom_point(color="black", pch=21) + 
    scale_fill_viridis(option="B") + 
    xlim(0, 1) +
    ylim(0, 1) +
    theme(legend.position="right") + 
    theme_minimal() + 
    geom_vline(xintercept=0.67)
  
  ggplot(subsampled, aes(CD27, HLA_DR, fill=timepoint)) + 
    geom_point(color="black", pch=21) + 
    xlim(0, 1) +
    ylim(0, 1) +
    theme(legend.position="right") + 
    theme_minimal() + 
    geom_vline(xintercept=0.67)
  
  n.sample <- 200000
  set.seed(888)
  subsampled <- dt[non_B==F, .SD[sample(.N, n.sample)], by=.(timepoint)] 
  subsampled$density <- get_density(subsampled$CD27, subsampled$CD19, n=80)
  
  ggplot(subsampled[meta=="CD11c_Effector"], aes(CD27, CD19)) + 
    geom_point(aes(color=density), stroke=0, size=3, alpha=0.8) + 
    scale_color_viridis(option='magma', name="Density", direction=-1) + 
    ylim(0.5,1) + 
    theme(legend.position="right") + 
    theme_minimal() + 
    geom_vline(xintercept=0.67)
  
  dt[, `:=`(CD27_pos=F)]
  
  dt[CD27>0.67, CD27_pos:=T]
  
  dt[, .N, by=.(CD27_pos)]
  
  n.sample <- 2000
  set.seed(888)
  subsampled <- dt[non_B==F, .SD[sample(.N, n.sample)], by=.(timepoint)]
  
  ggplot(subsampled, aes(CD27, HLA_DR, fill=timepoint)) + 
    geom_point(color="black", pch=21) + 
    xlim(0, 1.1) +
    ylim(0, 1.1) +
    theme(legend.position="right") + 
    theme_minimal() + 
    geom_vline(xintercept=0.67)
  ggsave(paste0(pa, "CD27_pos_B_cell_biaxial.eps"), width=8, height=6)
  ggsave(paste0(pa, "CD27_pos_B_cell_biaxial.png"), width=8, height=6)
  
  subsampled$density <- get_density(subsampled$CD27, subsampled$CD19, n=80)
  ggplot(subsampled, aes(CD27, CD19)) + 
    geom_point(aes(color=density), stroke=0, size=3, alpha=0.8) + 
    scale_color_viridis(option='magma', name="Density", direction=-1) + 
    xlim(0, 1.1) +
    ylim(0, 1.1) +
    theme(legend.position="right") + 
    theme_minimal() + 
    geom_vline(xintercept=0.67)
  ggsave(paste0(pa, "CD27_pos_density_B_cell_biaxial.eps"), width=8, height=6)
  ggsave(paste0(pa, "CD27_pos_density_B_cell_biaxial.png"), width=8, height=6)
  
  return(dt)
}


bMarkerLongitudinalPlots <- function(dt, 
                                     cc=subset.colors, 
                                     gc=group.colors, 
                                     tc=time.colors,
                                     pa=images.path) {
  # Generates box and violin plots of longitudinal marker expression and positive expression frequencies 
  # by subset and total B cells across flu year1 vaccine time points
  # Inputs:
  #   dt - data.table of B cells from PBMC flow
  #   cc - named vector of subset colors
  #   gc - named vector of group colors
  #   tc - named vector of time point colors
  #   pa - path to images folder
  # Outputs:
  #   png of marker expression violin plots and marker-positive cell frequency box plots
  
  #CXCR3 frequency
  n.sample <- 2000
  set.seed(888)
  subsampled <- dt[, .SD[sample(.N, n.sample)], by=.(timepoint)] 
  
  ggplot(subsampled, aes(CD183, HLA_DR, fill=CXCR3_pos)) + 
    geom_point(color="black", pch=21) + 
    xlim(0.3, 1) +
    theme(legend.position="right") + 
    theme_minimal() + 
    geom_vline(xintercept=0.85)
  ggsave(paste0(pa, "Perc_CXCR3_pos_biaxial_plot.png"), width=12, height=12)
  
  percxcr3stat <- table(dt[, .(donor, CXCR3_pos, timepoint, group)]) %>%
    as.data.table() %>%
    .[CXCR3_pos==T, NCXCR3:=N, by=.(donor, timepoint)] %>%
    .[, PCXCR3:=100*NCXCR3/sum(N), by=.(donor, timepoint)] %>%
    .[, c("N", "NCXCR3"):=NULL] %>%
    .[CXCR3_pos==T] %>%
    .[, CXCR3_pos:=NULL] %>% 
    .[PCXCR3!=0.000000] %>%
    .[, donor:=factor(donor, levels=donors)] %>%
    .[, group:=factor(group, levels=groups)] %>% 
    .[, timepoint:=factor(timepoint, levels=times)]
  
  ggplot(percxcr3stat, aes(timepoint, PCXCR3)) + 
    facet_grid(~group) + 
    geom_boxplot(aes(fill=group), width=0.4, color="black", position="dodge") + 
    geom_line(aes(group=donor), size=1) +
    geom_point(size=8, color='black') +
    ylab("Percent CXCR3+ B cells (%)") +
    xlab("Timepoint") +
    scale_fill_manual(values=gc) +
    ylim(0, 10) + 
    theme_minimal() +
    theme(text=element_text(size=35), 
          legend.position="none")
  ggsave(paste0(pa, "Perc_CXCR3_timepoint.png"), width=14, height=12)
  
  for (gr in groups) {
    
  percxcr3stat.sub <- percxcr3stat[group==gr]
  
  gt.vec <- percxcr3stat.sub[timepoint=="day0", PCXCR3]
  tc.vec <- percxcr3stat.sub[timepoint=="day7", PCXCR3]
  stat.test <- wilcox.test(gt.vec, tc.vec, paired=TRUE, alternative="less")
  
  ggplot(percxcr3stat.sub[timepoint!="day90"], aes(timepoint, PCXCR3))  + 
    geom_boxplot(aes(fill=group), width=0.4, color="black", position="dodge") + 
    geom_line(aes(group=donor), size=1) +
    geom_point(size=8, color='black') +
    scale_fill_manual(values=gc) +
    ggtitle(gr) +
    geom_signif(comparisons=list(c("day0", "day7")), 
                annotations=(paste("p=", round(stat.test$p.value, digits=2))), 
                map_signif_level=TRUE, size=0.6, textsize=12, tip_length=0, vjust=0.1) + 
    ylim(0, 10) + 
    ylab("Percent CXCR3+ B cells (%)") +
    xlab("") +
    theme_minimal() +
    theme(text=element_text(size=35), 
          legend.position='none')
  ggsave(paste0(pa, gr, "_Perc_CXCR3_timepoint_stats.eps"), width=10, height=12)
  ggsave(paste0(pa, gr, "_Perc_CXCR3_timepoint_stats.png"), width=10, height=12)
  
  }
  
  #CD69 frequency
  n.sample <- 2000
  set.seed(888)
  subsampled <- dt[, .SD[sample(.N, n.sample)], by=.(timepoint)] 
  
  ggplot(subsampled, aes(CD69, HLA_DR, fill=CD69_pos)) + 
    geom_point(color="black", pch=21) + 
    xlim(0.3, 1) +
    theme(legend.position="right") + 
    theme_minimal() + 
    geom_vline(xintercept=0.45)
  ggsave(paste0(pa, "Perc_CD69_pos_biaxial_plot.png"), width=12, height=12)
  
  perCD69stat <- table(dt[, .(donor, CD69_pos, timepoint, group)]) %>%
    as.data.table() %>%
    .[CD69_pos==T, NCD69:=N, by=.(donor, timepoint)] %>%
    .[, PCD69:=100*NCD69/sum(N), by=.(donor, timepoint)] %>%
    .[, c("N", "NCD69"):=NULL] %>%
    .[CD69_pos==T] %>%
    .[, CD69_pos:=NULL] %>% 
    .[PCD69!=0.000000] %>%
    .[, donor:=factor(donor, levels=donors)] %>%
    .[, group:=factor(group, levels=groups)] %>% 
    .[, timepoint:=factor(timepoint, levels=times)]
  
  ggplot(perCD69stat, aes(timepoint, PCD69)) + 
    facet_grid(~group) + 
    geom_boxplot(aes(fill=group), width=0.4, color="black", position="dodge") + 
    geom_line(aes(group=donor), size=1) +
    geom_point(size=8, color='black') +
    scale_fill_manual(values=gc) +
    ylab("Percent CD69+ B cells (%)") +
    xlab("Timepoint") +
    ylim(0, 75) + 
    theme_minimal() +
    theme(text=element_text(size=35), 
          legend.position="none")
  ggsave(paste0(pa, "Perc_CD69_timepoint.png"), width=14, height=12)
  
  for (gr in groups) {
    
    perCD69stat.sub <- perCD69stat[group==gr]
    
    gt.vec <- perCD69stat.sub[timepoint=="day0", PCD69]
    tc.vec <- perCD69stat.sub[timepoint=="day7", PCD69]
    stat.test <- wilcox.test(gt.vec, tc.vec, paired=TRUE, alternative="less")
    
    ggplot(perCD69stat.sub[timepoint!="day90"], aes(timepoint, PCD69))  + 
      geom_boxplot(aes(fill=group), width=0.4, color="black", position="dodge") + 
      geom_line(aes(group=donor), size=1) +
      geom_point(size=8, color='black') +
      scale_fill_manual(values=gc) +
      ggtitle(gr) +
      geom_signif(comparisons=list(c("day0", "day7")), 
                  annotations=(paste("p=", round(stat.test$p.value, digits=2))), 
                  map_signif_level=TRUE, size=0.6, textsize=12, tip_length=0, vjust=0.1) + 
      ylim(0, 75) + 
      ylab("Percent CD69+ B cells (%)") +
      xlab("") +
      theme_minimal() +
      theme(text=element_text(size=35), 
            legend.position='none')
    ggsave(paste0(pa, gr, "_Perc_CD69_timepoint_stats.eps"), width=10, height=12)
    ggsave(paste0(pa, gr, "_Perc_CD69_timepoint_stats.png"), width=10, height=12)
    
  }
  
  #activation marker expression levels
  
  markers <- c("CD19", "CD69", "CD183", "CD20", "CD11c")
  
  n.sample <- 3000
  set.seed(888)
  sub <- dt[, .SD[sample(.N, n.sample)], by=.(group)] %>% 
    as.data.table() %>%
    .[, timepoint:=factor(timepoint, levels=times)] %>%
    .[, group:=factor(group, levels=groups)]
  
  for (mar in markers) {
    
    ggplot(sub, aes(timepoint, eval(parse(text=mar)), color=timepoint)) + 
      facet_grid(~group) + 
      geom_quasirandom(alpha=1.0, size=3, width=0.3) +
      scale_color_manual(values=tc) +
      stat_summary(fun="median", aes(fill="Median"), geom="point", shape=18, size=6, 
                   position=position_dodge(), color='black', na.rm=T) + 
      theme_minimal() +
      scale_fill_manual(values="black", name="") + 
      ylab(paste(mar, "expression")) +
      xlab(" ") +
      ggtitle("Total B cells") + 
      theme(text=element_text(size=30), 
            axis.text.x=element_text(angle=90, hjust=1, size=20),
            legend.position='none')
    ggsave(paste0(pa, mar, "_expr_Bcells_time_bee.png"), height=8, width=12)
    
    med <- dt[timepoint!="day90", lapply(.SD, FUN=median), .SDcols=mar, by=.(group, timepoint)] %>%
      .[timepoint=="day0", median.0:=eval(parse(text=mar)), by=group] %>%
      .[timepoint=="day7", median.7:=eval(parse(text=mar)), by=group] %>%
      .[, timepoint:=NULL] %>%
      .[, (text=mar):=NULL] %>%
      melt(id.vars="group", variable.name="time", value.name= "median") %>% 
      as.data.table() %>% 
      .[median!="NA"]
    FC <- acast(med, group~time, value.var='median') %>% 
      as.data.table(keep.rownames="group") %>% 
      .[, group:=factor(group, levels=groups)] %>%
      .[, log2FC:=log2(median.7/median.0)]
    
    ggplot(FC, aes(group, log2FC, fill=group)) + 
      geom_col(width=0.4, color="black") + 
      theme_minimal() +
      scale_fill_manual(values=gc) +
      ylab(paste(mar, "log2 fold change day7/day0")) +
      xlab(" ") +
      ggtitle("Total B cells") + 
      theme(text=element_text(size=30), 
            axis.text.x=element_text(angle=90, hjust=1, size=20),
            legend.text=element_text(size=20),
            legend.position='none')
    ggsave(paste0(pa, mar, "_Bcells_day7_FC_group.png"), height=10, width=8)
    
    medians <- dt[timepoint!="day90", lapply(.SD, FUN=median), .SDcols=mar, by=.(group, timepoint, donor)] 
    
    ggplot(medians, aes(timepoint, eval(parse(text=mar)), color=timepoint, shape=group)) + 
      facet_grid(~group) + 
      geom_line(aes(group=donor, color='black'), size=0.5) +
      geom_point(size=8) + 
      stat_summary(fun="median", aes(fill="Median"), geom="crossbar", size=1, width=0.5, 
                   position=position_dodge(), color='black', na.rm=T) + 
      scale_color_manual(values=tc) +
      theme_minimal() +
      ylab(paste(mar, "median expression")) +
      xlab(" ") +
      ggtitle("Total B cells") + 
      theme(text=element_text(size=30), 
            axis.text.x=element_text(angle=90, hjust=1, size=20),
            legend.text=element_text(size=20),
            legend.position='none')
    ggsave(paste0(pa, mar, "_expr_Bcells_time_medians.png"), height=10, width=8)
    
  }
  
  #activation marker expression levels by B cell subset
  
  markers <- c("CD19", "CD69", "CD183", "CD20", "CD11c")
  
  n.sample <- 100000
  set.seed(888)
  sub <- dt[, .SD[sample(.N, n.sample)], by=.(group)] %>% 
    as.data.table() %>%
    .[, timepoint:=factor(timepoint, levels=times)] %>%
    .[, group:=factor(group, levels=groups)]
  
  for (mar in markers) {
    for (subs in subsets) {
    
    ggplot(sub[meta==subs], aes(timepoint, eval(parse(text=mar)), color=timepoint)) + 
      facet_grid(~group) + 
      geom_quasirandom(alpha=1.0, size=3, width=0.3) +
      scale_color_manual(values=tc) +
      stat_summary(fun="median", aes(fill="Median"), geom="point", shape=18, size=6, 
                   position=position_dodge(), color='black', na.rm=T) + 
      theme_minimal() +
      scale_fill_manual(values="black", name="") + 
      ylab(paste(mar, "expression")) +
      xlab(" ") +
      ggtitle(paste(subs, "\nB cells")) + 
      theme(text=element_text(size=30), 
            axis.text.x=element_text(angle=90, hjust=1, size=25),
            legend.position='none')
    ggsave(paste0(pa, mar, "_expr_", subs, "_Bcells_time_bee.png"), height=8, width=14)
    
    medians <- dt[meta==subs & timepoint!="day90", lapply(.SD, FUN=median), .SDcols=mar, by=.(group, timepoint, donor)] 
    
    ggplot(medians, aes(timepoint, eval(parse(text=mar)), color=timepoint, shape=group)) + 
      facet_grid(~group) + 
      geom_line(aes(group=donor, color='black'), size=0.5) +
      geom_point(size=8) + 
      stat_summary(fun="median", aes(fill="Median"), geom="crossbar", size=1, width=0.5, 
                   position=position_dodge(), color='black', na.rm=T) + 
      scale_color_manual(values=tc) +
      theme_minimal() +
      ylab(paste(mar, "median expression")) +
      xlab(" ") +
      ggtitle(paste(subs, "\nB cells")) + 
      theme(text=element_text(size=30), 
            axis.text.x=element_text(angle=90, hjust=1, size=25),
            legend.position='none')
    ggsave(paste0(pa, mar, "_expr_", subs, "_Bcells_time_medians.png"), height=10, width=8)
    
    }
  }
  
  #CD27- Effector B cell subset only
  
  markers <- c("CD19", "CD69", "CD183", "CD20", "CD11c", "HLA_DR")
 
  n.sample <- 4000
  set.seed(888)
  sub <- dt[meta=="CD27_neg_Effector" & timepoint!="day90", .SD[sample(.N, n.sample)], by=.(group)] %>% 
    as.data.table() %>%
    .[, timepoint:=factor(timepoint)] %>%
    .[, group:=factor(group, levels=groups)]
  
  for (mar in markers) {
    
    ggplot(sub, aes(timepoint, eval(parse(text=mar)), color=timepoint)) + 
      facet_grid(~group) + 
      geom_quasirandom(alpha=1.0, size=3, width=0.3) +
      scale_color_manual(values=tc) +
      stat_summary(fun="median", aes(fill="Median"), geom="point", shape=18, size=6, 
                   position=position_dodge(), color='black', na.rm=T) + 
      theme_minimal() +
      scale_fill_manual(values="black", name="") + 
      ylab(paste(mar, "expression")) +
      xlab(" ") +
      ggtitle("CD27- Effector B cells") + 
      theme(text=element_text(size=30), 
            axis.text.x=element_text(angle=90, hjust=1, size=30),
            legend.position='none')
    ggsave(paste0(pa, mar, "_expr_CD27neg_EffectorB_time_bee.png"), height=8, width=12)
 
    medians <- dt[meta=="CD27_neg_Effector" & timepoint!="day90", lapply(.SD, FUN=median), .SDcols=mar, by=.(group, timepoint, donor)] 
    
    gt.vec.0 <- medians[timepoint=="day0" & group=="BR1", eval(parse(text=mar))]
    tc.vec.0 <- medians[timepoint=="day0" & group=="BR2", eval(parse(text=mar))]
    stat.test.0 <- wilcox.test(gt.vec.0, tc.vec.0, paired=FALSE, alternative="greater")
    
    gt.vec.7 <- medians[timepoint=="day7" & group=="BR1", eval(parse(text=mar))]
    tc.vec.7 <- medians[timepoint=="day7" & group=="BR2", eval(parse(text=mar))]
    stat.test.7 <- wilcox.test(gt.vec.7, tc.vec.7, paired=FALSE, alternative="greater")
    
    ggplot(medians[timepoint=="day0"], aes(group, eval(parse(text=mar)), fill=group, shape=group)) + 
      geom_boxplot(width=0.4, position="dodge") + 
      geom_point(color="black", size=8, stroke=1) +
      geom_signif(comparisons=list(c("BR1", "BR2")), 
                  annotations=(paste("p=", round(stat.test.0$p.value, digits=3))), 
                  map_signif_level=TRUE, size=0.6, textsize=12, tip_length=0, vjust=0.1) + 
      scale_fill_manual(values=gc) +
      theme_minimal() +
      ylim(0,1.1) +
      ylab(paste(mar, "median expression")) +
      xlab(" ") +
      ggtitle("Day 0 \n CD27- Effector B cells") + 
      theme(text=element_text(size=40), 
            title=element_text(size=32),
            legend.position='none')
    ggsave(paste0(pa, mar, "_CD27neg_EffectorB_day0_medians.png"), height=10, width=8)
    
    ggplot(medians[timepoint=="day7"], aes(group, eval(parse(text=mar)), fill=group, shape=group)) + 
      geom_boxplot(width=0.4, position="dodge") + 
      geom_point(color="black", size=8, stroke=1) +
      geom_signif(comparisons=list(c("BR1", "BR2")), 
                  annotations=(paste("p=", round(stat.test.7$p.value, digits=3))), 
                  map_signif_level=TRUE, size=0.6, textsize=12, tip_length=0, vjust=0.1) + 
      scale_fill_manual(values=gc) +
      theme_minimal() +
      ylim(0,1.1) +
      ylab(paste(mar, "median expression")) +
      xlab(" ") +
      ggtitle("Day 7 \n CD27- Effector B cells") + 
      theme(text=element_text(size=40), 
            title=element_text(size=32),
            legend.position='none')
    ggsave(paste0(pa, mar, "_CD27neg_EffectorB_day7_medians.png"), height=10, width=8)
  }

}


isotypeBTime <- function(dt, 
                        ic=isotype.colors, 
                        pa=images.path) {
  # Generates stacked bar graph with percentages of B cell isotype and condition
  # Inputs:
  #   dt - data.table of data
  #   ic - named vector of isotype colors
  #   pa - path to images folder
  # Outputs:
  #   eps of stacked bar plot with freq of B cell isotypes by timepoint
  
  temp <- dt[, .N, by=.(isotype, group, timepoint)] %>% 
    .[, P:=N*100/sum(N), by=.(group, timepoint)] %>%
    .[, isotype:=factor(isotype, levels=isotypes)] %>% 
    .[, timepoint:=factor(timepoint, levels=times)]
  
  ggplot(temp, aes(timepoint, P, fill=isotype)) + 
    geom_col(width=0.9, color="black", position="fill") + 
    facet_grid(~group) + 
    theme_minimal() +
    ylab("Percent of B cells (%)") +
    xlab(" ") +
    labs(fill="Isotype") +
    scale_fill_manual(values=ic) +
    theme(text=element_text(size=25), 
          axis.text.x=element_text(angle=270, size=20), 
          legend.position='right')
  ggsave(paste0(pa, "Stack_b_isotype.eps"), width=8, height=8)
  ggsave(paste0(pa, "Stack_b_isotype.png"), width=8, height=8)

}


percentBoxBcellIso <- function(dt, 
                               ic=iso.colors, 
                               dc=donor.colors, 
                               tc=time.colors,
                               pa=images.path) {
  # Generates bar graph of B cell isotype percentages by timepoint and group
  # Inputs:
  #   dt - data.table  of B cells from PBMC flow
  #   ic - named vector of isotype colors
  #   dc - named vector of donor colors
  #   pa - path to images folder
  # Outputs:
  #   png of box plots with percent B cell isotypes by group and by subset + group
  
  for (iso in isotypes) {
    
    temp <- table(dt[, .(timepoint, donor, isotype, group)]) %>%
      as.data.table() %>%
      .[isotype==iso, NBcell:=N, by=.(timepoint, donor)] %>%
      .[, PBcell:=100*NBcell/sum(N), by=.(timepoint, donor)] %>%
      .[, c("N", "NBcell"):=NULL] %>% 
      .[isotype==iso] %>%
      .[PBcell!="NaN"] %>%
      .[PBcell!=0.00000] %>% 
      .[, timepoint:=factor(timepoint, levels=times)]
    
    ggplot(temp[timepoint!="day90"], aes(timepoint, PBcell)) + 
      geom_boxplot(aes(fill=timepoint), width=0.4, position="dodge") + 
      geom_point(color="black", size=6, stroke=1) +
      scale_fill_manual(values=tc, name="Timepoint") + 
      facet_wrap(~group) +
      ylim(0,100) +
      ylab("Percent of B cells (%)") +
      xlab("") +
      ggtitle(paste(iso)) + 
      theme_minimal() +
      theme(text=element_text(size=30),
            legend.position='none')
    ggsave(paste0(pa, iso, "_perc_Bcell_donor_day_boxplot.png"), width=12, height=8)
    
    temp.BR1 <- temp[group=="BR1"]
    
    gt.vec <- temp.BR1[timepoint=="day0", PBcell]
    tc.vec <- temp.BR1[timepoint=="day7", PBcell]
    stat.test.iso.BR1 <- wilcox.test(gt.vec, tc.vec, paired=TRUE, alternative="l")
    
    ggplot(temp.BR1[timepoint!="day90"], aes(timepoint, PBcell))  + 
      geom_boxplot(aes(fill=timepoint), width=0.4, position="dodge") + 
      geom_point(color="black", size=6, stroke=1) +
      scale_fill_manual(values=tc, name="Timepoint") + 
      geom_signif(comparisons=list(c("day0", "day7")), 
                  annotations=(paste("p=", round(stat.test.iso.BR1$p.value, digits=2))), 
                  map_signif_level=TRUE, size=0.6, textsize=10, tip_length=0, vjust=0.1) + 
      ylim(0,100) +
      ylab(paste("Percent", iso, "of B cells (%)")) +
      xlab("") +
      ggtitle(paste("BR1 -", iso)) + 
      theme_minimal() +
      theme(text=element_text(size=20), 
            legend.position='none')
    ggsave(paste0(pa, iso, "_BR1_isotype_perc_day_stats.png"), width=6, height=8)
    
    temp.BR2 <- temp[group=="BR2"]
    
    gt.vec <- temp.BR2[timepoint=="day0", PBcell]
    tc.vec <- temp.BR2[timepoint=="day7", PBcell]
    stat.test.iso.BR2 <- wilcox.test(gt.vec, tc.vec, paired=TRUE, alternative="l")
    
    ggplot(temp.BR2[timepoint!="day90"], aes(timepoint, PBcell))  + 
      geom_boxplot(aes(fill=timepoint), width=0.4, position="dodge") + 
      geom_point(color="black", size=6, stroke=1) +
      scale_fill_manual(values=tc, name="Timepoint") + 
      geom_signif(comparisons=list(c("day0", "day7")), 
                  annotations=(paste("p=", round(stat.test.iso.BR2$p.value, digits=2))), 
                  map_signif_level=TRUE, size=0.6, textsize=10, tip_length=0, vjust=0.1) + 
      ylim(0,100) +
      ylab(paste("Percent", iso, "of B cells (%)")) +
      xlab("") +
      ggtitle(paste("BR2 -", iso)) + 
      theme_minimal() +
      theme(text=element_text(size=20), 
            legend.position='none')
    ggsave(paste0(pa, iso, "_BR2_isotype_perc_day_stats.png"), width=6, height=8)
  }
  
  
  #CD27- Effector B cells
  
  for (iso in c("IgG", "IgA")) {
    
    temp <- table(dt[meta=="CD27_neg_Effector" & timepoint!="day90", .(timepoint, donor, isotype, group)]) %>%
      as.data.table() %>%
      .[isotype==iso, NBcell:=N, by=.(timepoint, donor)] %>%
      .[, PBcell:=100*NBcell/sum(N), by=.(timepoint, donor)] %>%
      .[, c("N", "NBcell"):=NULL] %>% 
      .[isotype==iso] %>%
      .[PBcell!="NaN"] %>%
      .[PBcell!=0.00000000] %>% 
      .[, timepoint:=factor(timepoint, levels=times)] %>%
      .[, group:=factor(group, levels=groups)]
    
    temp.day0 <- temp[timepoint=="day0"]
    
    gt.vec <- temp.day0[group=="BR1", PBcell]
    tc.vec <- temp.day0[group=="BR2", PBcell]
    stat.test.iso.day0 <- wilcox.test(gt.vec, tc.vec, paired=F, alternative="g")
    
    ggplot(temp.day0[timepoint!="day90"], aes(group, PBcell))  + 
      geom_boxplot(aes(fill=group), width=0.3, color="black", position="dodge") + 
      geom_point(size=6, color="black") +
      geom_signif(comparisons=list(c("BR1", "BR2")), 
                  annotations=(paste("p=", round(stat.test.iso.day0$p.value, digits=2))), 
                  map_signif_level=TRUE, size=0.6, textsize=10, tip_length=0, vjust=0.1) + 
      ylim(0,100) + 
      ylab(paste("Percentage", iso, "cells (%)")) +
      xlab("") +
      ggtitle(paste0("Day 0 \n", iso, "+ \nCD27- Effector B cells")) + 
      scale_fill_manual(values=group.colors, name="") +
      theme_minimal() +
      theme(text=element_text(size=32), 
            title=element_text(size=28), 
            legend.position='none')
    ggsave(paste0(pa, iso, "_CD27neg_EffectorB_day0_perc_stats.png"), width=6, height=8)
    
    temp.day7 <- temp[timepoint=="day7"]
    
    gt.vec <- temp.day7[group=="BR1", PBcell]
    tc.vec <- temp.day7[group=="BR2", PBcell]
    stat.test.iso.day7 <- wilcox.test(gt.vec, tc.vec, paired=F, alternative="g")
    
    ggplot(temp.day7[timepoint!="day90"], aes(group, PBcell))  + 
      geom_boxplot(aes(fill=group), width=0.3, color="black", position="dodge") + 
      geom_point(size=6, color="black") +
      geom_signif(comparisons=list(c("BR1", "BR2")), 
                  annotations=(paste("p=", round(stat.test.iso.day7$p.value, digits=2))), 
                  map_signif_level=TRUE, size=0.6, textsize=10, tip_length=0, vjust=0.1) +
      ylim(0,100) + 
      ylab(paste("Percentage", iso, "cells (%)")) +
      xlab("") +
      ggtitle(paste0("Day 7 \n", iso, "+ \nCD27- Effector B cells")) + 
      scale_fill_manual(values=group.colors, name="") +
      theme_minimal() +
      theme(text=element_text(size=32), 
            title=element_text(size=28), 
            legend.position='none')
    ggsave(paste0(pa, iso, "_CD27neg_EffectorB_day7_perc_stats.png"), width=6, height=8)
  }
  
  for (iso in c("IgMD")) {
    
    temp <- table(dt[meta=="CD27_neg_Effector" & timepoint!="day90", .(timepoint, donor, isotype, group)]) %>%
      as.data.table() %>%
      .[isotype==iso, NBcell:=N, by=.(timepoint, donor)] %>%
      .[, PBcell:=100*NBcell/sum(N), by=.(timepoint, donor)] %>%
      .[, c("N", "NBcell"):=NULL] %>% 
      .[isotype==iso] %>%
      .[PBcell!="NaN"] %>%
      .[PBcell!=0.00000000] %>% 
      .[, timepoint:=factor(timepoint, levels=times)] %>%
      .[, group:=factor(group, levels=groups)]
    
    temp.day0 <- temp[timepoint=="day0"]
    
    gt.vec <- temp.day0[group=="BR1", PBcell]
    tc.vec <- temp.day0[group=="BR2", PBcell]
    stat.test.iso.day0 <- wilcox.test(gt.vec, tc.vec, paired=F, alternative="l")
    
    ggplot(temp.day0[timepoint!="day90"], aes(group, PBcell))  + 
      geom_boxplot(aes(fill=group), width=0.3, color="black", position="dodge") + 
      geom_point(size=6, color="black") +
      geom_signif(comparisons=list(c("BR1", "BR2")), 
                  annotations=(paste("p=", round(stat.test.iso.day0$p.value, digits=2))), 
                  map_signif_level=TRUE, size=0.6, textsize=10, tip_length=0, vjust=0.1) + 
      ylim(0,100) + 
      ylab(paste("Percentage", iso, "cells (%)")) +
      xlab("") +
      ggtitle(paste0("Day 0 \n", iso, "+ \nCD27- Effector B cells")) + 
      scale_fill_manual(values=group.colors, name="") +
      theme_minimal() +
      theme(text=element_text(size=32), 
            title=element_text(size=28), 
            legend.position='none')
    ggsave(paste0(pa, iso, "_CD27neg_EffectorB_day0_perc_stats.png"), width=6, height=8)
    
    temp.day7 <- temp[timepoint=="day7"]
    
    gt.vec <- temp.day7[group=="BR1", PBcell]
    tc.vec <- temp.day7[group=="BR2", PBcell]
    stat.test.iso.day7 <- wilcox.test(gt.vec, tc.vec, paired=F, alternative="l")
    
    ggplot(temp.day7[timepoint!="day90"], aes(group, PBcell))  + 
      geom_boxplot(aes(fill=group), width=0.3, color="black", position="dodge") + 
      geom_point(size=6, color="black") +
      geom_signif(comparisons=list(c("BR1", "BR2")), 
                  annotations=(paste("p=", round(stat.test.iso.day7$p.value, digits=2))), 
                  map_signif_level=TRUE, size=0.6, textsize=10, tip_length=0, vjust=0.1) +
      ylim(0,100) + 
      ylab(paste("Percentage", iso, "cells (%)")) +
      xlab("") +
      ggtitle(paste0("Day 7 \n", iso, "+ \nCD27- Effector B cells")) + 
      scale_fill_manual(values=group.colors, name="") +
      theme_minimal() +
      theme(text=element_text(size=32),               
            title=element_text(size=28), 
            legend.position='none')
    ggsave(paste0(pa, iso, "_CD27neg_EffectorB_day7_perc_stats.png"), width=6, height=8)
  }
  
  
  # by B cell subset
  
  for (subs in c("CD27_neg_Effector", "Core_Memory")) {
    for (iso in c("IgG", "IgA", "IgMD")) {
      
      temp <- table(dt[meta==subs, .(timepoint, donor, isotype, group)]) %>%
        as.data.table() %>%
        .[isotype==iso, NBcell:=N, by=.(timepoint, donor)] %>%
        .[, PBcell:=100*NBcell/sum(N), by=.(timepoint, donor)] %>%
        .[, c("N", "NBcell"):=NULL] %>% 
        .[isotype==iso] %>%
        .[PBcell!="NaN"] %>%
        .[PBcell!=0.00000000] %>% 
        .[, timepoint:=factor(timepoint, levels=times)]
      
      ggplot(temp[timepoint!="day90"], aes(timepoint, PBcell)) + 
        geom_boxplot(aes(fill=timepoint), width=0.4, position="dodge") + 
        geom_point(color="black", size=6, stroke=1) +
        scale_fill_manual(values=tc, name="Timepoint") + 
        facet_wrap(~group) +
        ylim(0,100) +
        ylab(paste("Percent of", subs, "cells (%)")) +
        xlab("") +
        ggtitle(paste(iso, subs, "cells")) + 
        theme_minimal() +
        theme(text=element_text(size=20), 
              legend.position='none')
      ggsave(paste0(pa, iso, "_", subs, "_perc_Bcell_donor_day_boxplot.png"), width=12, height=8)
      
      temp.BR1 <- temp[group=="BR1"]
      
      gt.vec <- temp.BR1[timepoint=="day0", PBcell]
      tc.vec <- temp.BR1[timepoint=="day7", PBcell]
      stat.test.iso.BR1 <- wilcox.test(gt.vec, tc.vec, paired=TRUE, alternative="l")
      
      ggplot(temp.BR1[timepoint!="day90"], aes(timepoint, PBcell))  + 
        geom_boxplot(aes(fill=timepoint), width=0.4, position="dodge") + 
        geom_point(color="black", size=6, stroke=1) +
        scale_fill_manual(values=tc, name="Timepoint") + 
        geom_signif(comparisons=list(c("day0", "day7")), 
                    annotations=(paste("p=", round(stat.test.iso.BR1$p.value, digits=3))), 
                    map_signif_level=TRUE, size=0.6, textsize=10, tip_length=0, vjust=0.1) + 
        ylim(0,100) +
        ylab(paste("Percent", iso, "of", subs, "cells (%)")) +
        xlab("") +
        ggtitle(paste("BR1 -", iso, subs)) + 
        theme_minimal() +
        theme(text=element_text(size=16), 
              legend.position='none')
      ggsave(paste0(pa, iso, "_", subs, "_BR1_isotype_perc_day_stats.png"), width=6, height=8)
      
      temp.BR2 <- temp[group=="BR2"]
      
      gt.vec <- temp.BR2[timepoint=="day0", PBcell]
      tc.vec <- temp.BR2[timepoint=="day7", PBcell]
      stat.test.iso.BR2 <- wilcox.test(gt.vec, tc.vec, paired=TRUE, alternative="l")
      
      ggplot(temp.BR2[timepoint!="day90"], aes(timepoint, PBcell))  + 
        geom_boxplot(aes(fill=timepoint), width=0.4, position="dodge") + 
        geom_point(color="black", size=6, stroke=1) +
        scale_fill_manual(values=tc, name="Timepoint") + 
        geom_signif(comparisons=list(c("day0", "day7")), 
                    annotations=(paste("p=", round(stat.test.iso.BR2$p.value, digits=3))), 
                    map_signif_level=TRUE, size=0.6, textsize=10, tip_length=0, vjust=0.1) + 
        ylim(0,100) +
        ylab(paste("Percent", iso, "of", subs, "cells (%)")) +
        xlab("") +
        ggtitle(paste("BR2 -", iso, subs)) + 
        theme_minimal() +
        theme(text=element_text(size=16), 
              legend.position='none')
      ggsave(paste0(pa, iso, "_", subs, "_BR2_isotype_perc_day_stats.png"), width=6, height=8)
    }
  }
  
  for (subs in c("CD95_Memory", "Plasma")) {
    for (iso in c("IgG")) {
      
      temp <- table(dt[meta==subs, .(timepoint, donor, isotype, group)]) %>%
        as.data.table() %>%
        .[isotype==iso, NBcell:=N, by=.(timepoint, donor)] %>%
        .[, PBcell:=100*NBcell/sum(N), by=.(timepoint, donor)] %>%
        .[, c("N", "NBcell"):=NULL] %>% 
        .[isotype==iso] %>%
        .[PBcell!="NaN"] %>%
        .[PBcell!=0.00000000] %>% 
        .[, timepoint:=factor(timepoint, levels=times)]
      
      ggplot(temp[timepoint!="day90"], aes(timepoint, PBcell)) + 
        geom_boxplot(aes(fill=timepoint), width=0.4, position="dodge") + 
        geom_point(color="black", size=6, stroke=1) + 
        facet_wrap(~group) +
        ylim(0,100) +
        ylab(paste("Percent of", subs, "cells (%)")) +
        xlab("") +
        ggtitle(paste(iso, subs, "cells")) + 
        scale_fill_manual(values=tc, name="Timepoint") + 
        theme_minimal() +
        theme(text=element_text(size=20), 
              legend.position="none")
      ggsave(paste0(pa, iso, "_", subs, "_perc_Bcell_donor_day_boxplot.png"), width=10, height=8)
      
      temp.BR1 <- temp[group=="BR1"]
      
      gt.vec <- temp.BR1[timepoint=="day0", PBcell]
      tc.vec <- temp.BR1[timepoint=="day7", PBcell]
      stat.test.iso.BR1 <- wilcox.test(gt.vec, tc.vec, paired=TRUE, alternative="l")
      
      ggplot(temp.BR1[timepoint!="day90"], aes(timepoint, PBcell))  + 
        geom_boxplot(aes(fill=timepoint), width=0.4, position="dodge") + 
        geom_point(color="black", size=6, stroke=1) +
        geom_signif(comparisons=list(c("day0", "day7")), 
                    annotations=(paste("p=", round(stat.test.iso.BR1$p.value, digits=3))), 
                    map_signif_level=TRUE, size=0.6, textsize=10, tip_length=0, vjust=0.1) + 
        ylim(0,100) +
        ylab(paste("Percent", iso, "of", subs, "cells (%)")) +
        xlab("") +
        ggtitle(paste("BR1 -", iso, subs)) + 
        scale_fill_manual(values=tc, name="Timepoint") + 
        theme_minimal() +
        theme(text=element_text(size=20), 
              legend.position='none')
      ggsave(paste0(pa, iso, "_", subs, "_BR1_isotype_perc_day_stats.png"), width=6, height=10)
      
      temp.BR2 <- temp[group=="BR2"]
      
      gt.vec <- temp.BR2[timepoint=="day0", PBcell]
      tc.vec <- temp.BR2[timepoint=="day7", PBcell]
      stat.test.iso.BR2 <- wilcox.test(gt.vec, tc.vec, paired=TRUE, alternative="l")
      
      ggplot(temp.BR2[timepoint!="day90"], aes(timepoint, PBcell))  + 
        geom_boxplot(aes(fill=timepoint), width=0.4, position="dodge") + 
        geom_point(color="black", size=6, stroke=1) +
        scale_fill_manual(values=tc, name="Timepoint") + 
        geom_signif(comparisons=list(c("day0", "day7")), 
                    annotations=(paste("p=", round(stat.test.iso.BR2$p.value, digits=3))), 
                    map_signif_level=TRUE, size=0.6, textsize=10, tip_length=0, vjust=0.1) + 
        ylim(0,100) +
        ylab(paste("Percent", iso, "of", subs, "cells (%)")) +
        xlab("") +
        ggtitle(paste("BR2 -", iso, subs)) + 
        theme_minimal() +
        theme(text=element_text(size=20), 
              legend.position='none')
      ggsave(paste0(pa, iso, "_", subs, "_BR2_isotype_perc_day_stats.png"), width=6, height=10)
    }
  }

}


percentBoxBcellSub <- function(dt, 
                              cc=subset.colors, 
                              dc=donor.colors, 
                              tc=time.colors,
                              pa=images.path) {
  # Generates bar graph of B cell percentages by timepoint
  # Inputs:
  #   dt - data.table  of B cells from PBMC flow
  #   cc - named vector of subset colors
  #   dc - named vector of donor colors
  #   pa - path to images folder
  # Outputs:
  #   png of box plots with percent B cells by subset and donor
  
  dt <- dt[is.na(meta), meta:="non_B"] %>% 
    .[, `:=`(donor=factor(donor, levels=donors),
             cell_type=factor(cell_type), 
             timepoint=factor(timepoint, levels=times), 
             group=factor(group, levels=groups),
             isotype=factor(isotype, levels=isotypes), 
             meta=factor(meta, levels=c("CD27_pos_Effector", "CD27_neg_Effector",
                                        "Early_Memory", "Core_Memory", 
                                        "CD95_Memory", "T2MBC",
                                        "Transitional", "Naive", "Plasma", 
                                        "non_B")))]
  
  temp <- dt[meta!="non_B", .N, by=.(meta, group, timepoint)] %>% 
    .[, P:=N*100/sum(N), by=.(group, timepoint)] %>%
    .[, meta:=factor(meta, levels=subsets)] %>% 
    .[, timepoint:=factor(timepoint, levels=times)]
  
  ggplot(temp[timepoint!="day90"], aes(timepoint, P, fill=meta)) + 
    geom_col(width=0.8, color="black", position="fill") + 
    facet_grid(~group) + 
    theme_minimal() +
    ylab("Percent of B cells (%)") +
    xlab(" ") +
    labs(fill="Subset") +
    scale_fill_manual(values=cc) +
    theme(text=element_text(size=25), 
          axis.text.x=element_text(angle=270, size=20), 
          legend.position='right')
  ggsave(paste0(pa, "Stack_b_subset.eps"), width=10, height=10)
  ggsave(paste0(pa, "Stack_b_subset.png"), width=10, height=10)
  
  for (subs in subsets) {
    
    temp <- table(dt[, .(timepoint, donor, meta, group)]) %>%
      as.data.table() %>%
      .[meta==subs, NBcell:=N, by=.(timepoint, donor)] %>%
      .[, PBcell:=100*NBcell/sum(N), by=.(timepoint, donor)] %>%
      .[, c("N", "NBcell"):=NULL] %>% 
      .[meta==subs] %>%
      .[PBcell!="NaN"] %>%
      .[PBcell!=0.00000] %>% 
      .[, timepoint:=factor(timepoint, levels=times)]
    
    ggplot(temp[timepoint!="day90"], aes(timepoint, PBcell)) + 
      geom_boxplot(width=0.5, color="black", position="dodge") + 
      geom_point(aes(color=group), size=4, stroke=1) + 
      facet_wrap(~group) +
      ylab("Percent of live cells (%)") +
      xlab("") +
      ggtitle(paste(subs)) + 
      stat_summary(fun.y=mean, aes(fill="Mean"), geom="point", shape=21, size=4, stroke=1, 
                   color="black", position=position_dodge()) +
      scale_color_manual(values=group.colors, name="") +
      scale_fill_manual(values="light yellow", name="") + 
      theme_minimal() +
      theme(text=element_text(size=30), 
            axis.text.x=element_text(size=20, angle=270), 
            axis.ticks.x = element_line(size=25))
    ggsave(paste0(pa, subs, "_perc_Bcell_donor_day_boxplot.png"), width=12, height=8)
    
    temp.BR1 <- temp[group=="BR1"]

    gt.vec <- temp.BR1[timepoint=="day0", PBcell]
    tc.vec <- temp.BR1[timepoint=="day7", PBcell]
    stat.test.subs.BR1 <- wilcox.test(gt.vec, tc.vec, paired=TRUE, alternative="l")
    
    ggplot(temp.BR1[timepoint!="day90"], aes(timepoint, PBcell))  + 
      geom_boxplot(width=0.4, color="black", position="dodge") + 
      geom_point(size=4, aes(color=timepoint)) +
      #stat_summary(geom="crossbar", fun="mean", width=0.5) + 
      geom_signif(comparisons=list(c("day0", "day7")), 
                  annotations=(paste("p=", round(stat.test.subs.BR1$p.value, digits=2))), 
                  map_signif_level=TRUE, size=0.6, textsize=10, tip_length=0, vjust=0.1) + 
      ylab(paste("Percent", subs, "of live cells (%)")) +
      xlab("") +
      ggtitle(paste("BR1 -", subs)) + 
      scale_color_manual(values=tc, name="") +
      theme_minimal() +
      theme(text=element_text(size=20), 
            axis.text.x=element_text(hjust=1, angle=25),
            legend.text=element_text(size=20),
            legend.position='none')
    ggsave(paste0(pa, subs, "_BR1_subset_perc_day_stats.png"), width=6, height=8)
    
    temp.BR2 <- temp[group=="BR2"]
    
    gt.vec <- temp.BR2[timepoint=="day0", PBcell]
    tc.vec <- temp.BR2[timepoint=="day7", PBcell]
    stat.test.subs.BR2 <- wilcox.test(gt.vec, tc.vec, paired=TRUE, alternative="l")
    
    ggplot(temp.BR2[timepoint!="day90"], aes(timepoint, PBcell))  + 
      geom_boxplot(width=0.4, color="black", position="dodge") + 
      geom_point(size=4, aes(color=timepoint)) +
      #stat_summary(geom="crossbar", fun="mean", width=0.5) + 
      geom_signif(comparisons=list(c("day0", "day7")), 
                  annotations=(paste("p=", round(stat.test.subs.BR2$p.value, digits=2))), 
                  map_signif_level=TRUE, size=0.6, textsize=10, tip_length=0, vjust=0.1) + 
      ylab(paste("Percent", subs, "of live cells (%)")) +
      xlab("") +
      ggtitle(paste("BR2 -", subs)) + 
      scale_color_manual(values=tc, name="") +
      theme_minimal() +
      theme(text=element_text(size=20), 
            axis.text.x=element_text(hjust=1, angle=25),
            legend.text=element_text(size=20),
            legend.position='none')
    ggsave(paste0(pa, subs, "_BR2_subset_perc_day_stats.png"), width=6, height=8)
  }

}


makeBcellUmapSub <- function(dt=dat[non_B==F], 
                              mc=subset.colors, 
                              cc=time.colors, 
                              dc=donor.colors, 
                              sc=group.colors, 
                              ic=isotype.colors,
                              yc=active.colors,
                              pa=images.path,
                              fa=factors, 
                              su=surface.markers, 
                              markers=setdiff(colnames(dat), c(factors, active.factors, isotype.channels, "CD27_pos", "Viability", "Lin", "CD10", "CD23"))) {
  # Generates umaps of B cells
  # Inputs:
  #   dt - data.table with B cells
  #   cc - named vector of timepoint colors
  #   dc - named vector of donor colors
  #   yc - named vector of activation molecule positivity colors
  #   mc - named vector of b cell subset
  #   pa - path to images folder
  #   fa - vector of factor column names
  #   su - vector of surface marker column names
  #   markers - vector of marker column names 
  # Outputs:
  #   eps of surface marker-driven umap plots of B cells from PBMC data
  

    # subsample data by timepoint
    n.sample <- 12000
    set.seed(888)
    subsampled <- dt[, .SD[sample(.N, n.sample)], by=.(timepoint)]
    
    # setup and run umap on SUBSAMPLED data
    set.seed(888)
    umap.out <- umap(subsampled[, markers, with=F], n_neighbors=25, min_dist=0.3)
    subsampled[, `:=`(umap.1=umap.out[,1], umap.2=umap.out[,2])]
    
    # meta
    ggplot(subsampled, aes(umap.1, umap.2)) + 
      geom_density2d(color="#CCCCCC") +
      geom_point(size=2, aes(color=meta)) +
      theme_bw() + 
      labs(x=NULL, y=NULL) +
      xlim(min(subsampled$umap.1)-1, max(subsampled$umap.1)+1) + 
      ylim(min(subsampled$umap.2)-1.5, max(subsampled$umap.2)+1.5) +
      scale_color_manual(values=mc, name="Metaclusters") + 
      labs(x=NULL, y=NULL) +
      theme(panel.grid.major=element_blank(), 
            panel.grid.minor=element_blank(), 
            panel.border=element_blank(),
            axis.ticks=element_blank(), 
            axis.text=element_blank(), 
            legend.position="none")
    ggsave(paste0(pa, "umap/meta.png"))
    ggplot(subsampled, aes(umap.1, umap.2)) + 
      geom_density2d(color="#CCCCCC") +
      geom_point(size=2, aes(color=meta)) +
      theme_bw() + 
      labs(x=NULL, y=NULL) +
      xlim(min(subsampled$umap.1)-1, max(subsampled$umap.1)+1) + 
      ylim(min(subsampled$umap.2)-1.5, max(subsampled$umap.2)+1.5) +
      scale_color_manual(values=mc, name="Metaclusters") + 
      labs(x=NULL, y=NULL) +
      theme(panel.grid.major=element_blank(), 
            panel.grid.minor=element_blank(), 
            panel.border=element_blank(),
            legend.text=element_text(size=15), 
            axis.ticks=element_blank(), 
            axis.text=element_blank()) + 
      guides(color=guide_legend(override.aes = list(size=7, linetype=0)))
    ggsave(paste0(pa, "umap/meta_legend.png"))
    
    # donor
    ggplot(subsampled, aes(umap.1, umap.2)) + 
      geom_density2d(color="#CCCCCC") +
      geom_point(size=2, aes(color=donor)) +
      theme_bw() + 
      labs(x=NULL, y=NULL) +
      xlim(min(subsampled$umap.1)-1, max(subsampled$umap.1)+1) + 
      ylim(min(subsampled$umap.2)-1.5, max(subsampled$umap.2)+1.5) +
      scale_color_manual(values=dc, name="Subject") +
      labs(x=NULL, y=NULL) +
      theme(panel.grid.major=element_blank(), 
            panel.grid.minor=element_blank(), 
            panel.border=element_blank(),
            axis.ticks=element_blank(), 
            axis.text=element_blank(), 
            legend.position="none")
    ggsave(paste0(pa, "umap/subject.png"))
    ggplot(subsampled, aes(umap.1, umap.2)) + 
      geom_density2d(color="#CCCCCC") +
      geom_point(size=2, aes(color=donor)) +
      theme_bw() + 
      labs(x=NULL, y=NULL) +
      xlim(min(subsampled$umap.1)-1, max(subsampled$umap.1)+1) + 
      ylim(min(subsampled$umap.2)-1.5, max(subsampled$umap.2)+1.5) +
      scale_color_manual(values=dc, name="Subject") +
      labs(x=NULL, y=NULL) +
      theme(panel.grid.major=element_blank(), 
            panel.grid.minor=element_blank(), 
            panel.border=element_blank(),
            legend.text=element_text(size=15), 
            axis.ticks=element_blank(), 
            axis.text=element_blank())  + 
      guides(color=guide_legend(override.aes = list(size=5, linetype=0)))
    ggsave(paste0(pa, "umap/subject_legend.png"))
    
    # timepoint
    ggplot(subsampled, aes(umap.1, umap.2)) + 
      geom_density2d(color="#CCCCCC") +
      geom_point(size=2, aes(color=timepoint)) +
      theme_bw() + 
      labs(x=NULL, y=NULL) +
      xlim(min(subsampled$umap.1)-1, max(subsampled$umap.1)+1) + 
      ylim(min(subsampled$umap.2)-1.5, max(subsampled$umap.2)+1.5) +
      scale_color_manual(values=cc, name="timepoint") + 
      labs(x=NULL, y=NULL) +
      theme(panel.grid.major=element_blank(), 
            panel.grid.minor=element_blank(), 
            panel.border=element_blank(),
            axis.ticks=element_blank(), 
            axis.text=element_blank(), 
            legend.position="none")
    ggsave(paste0(pa, "umap/timepoint.png"))
    ggplot(subsampled, aes(umap.1, umap.2)) + 
      geom_density2d(color="#CCCCCC") +
      geom_point(size=2, aes(color=timepoint)) +
      theme_bw() + 
      labs(x=NULL, y=NULL) +
      xlim(min(subsampled$umap.1)-1, max(subsampled$umap.1)+1) + 
      ylim(min(subsampled$umap.2)-1.5, max(subsampled$umap.2)+1.5) +
      scale_color_manual(values=cc, name="timepoint") + 
      labs(x=NULL, y=NULL) +
      theme(panel.grid.major=element_blank(), 
            panel.grid.minor=element_blank(), 
            panel.border=element_blank(),
            legend.text=element_text(size=15), 
            axis.ticks=element_blank(), 
            axis.text=element_blank()) + 
      guides(color=guide_legend(override.aes = list(size=7, linetype=0)))
    ggsave(paste0(pa, "umap/timepoint_legend.png"))
    
   
    # group
    ggplot(subsampled, aes(umap.1, umap.2)) + 
      geom_density2d(color="#CCCCCC") +
      geom_point(size=2, aes(color=group)) +
      theme_bw() + 
      labs(x=NULL, y=NULL) +
      xlim(min(subsampled$umap.1)-1, max(subsampled$umap.1)+1) + 
      ylim(min(subsampled$umap.2)-1.5, max(subsampled$umap.2)+1.5) +
      scale_color_manual(values=sc, name="group") + 
      labs(x=NULL, y=NULL) +
      theme(panel.grid.major=element_blank(), 
            panel.grid.minor=element_blank(), 
            panel.border=element_blank(),
            axis.ticks=element_blank(), 
            axis.text=element_blank(), 
            legend.position="none")
    ggsave(paste0(pa, "umap/group.png"))
    ggplot(subsampled, aes(umap.1, umap.2)) + 
      geom_density2d(color="#CCCCCC") +
      geom_point(size=2, aes(color=group)) +
      theme_bw() + 
      labs(x=NULL, y=NULL) +
      xlim(min(subsampled$umap.1)-1, max(subsampled$umap.1)+1) + 
      ylim(min(subsampled$umap.2)-1.5, max(subsampled$umap.2)+1.5) +
      scale_color_manual(values=sc, name="group") + 
      labs(x=NULL, y=NULL) +
      theme(panel.grid.major=element_blank(), 
            panel.grid.minor=element_blank(), 
            panel.border=element_blank(),
            legend.text=element_text(size=15), 
            axis.ticks=element_blank(), 
            axis.text=element_blank()) + 
      guides(color=guide_legend(override.aes = list(size=7, linetype=0)))
    ggsave(paste0(pa, "umap/group_legend.png"))
    
  
    # isotype
    ggplot(subsampled, aes(umap.1, umap.2)) + 
      geom_density2d(color="#CCCCCC") +
      geom_point(size=2, aes(color=isotype)) +
      theme_bw() + 
      labs(x=NULL, y=NULL) +
      xlim(min(subsampled$umap.1)-1, max(subsampled$umap.1)+1) + 
      ylim(min(subsampled$umap.2)-1.5, max(subsampled$umap.2)+1.5) +
      scale_color_manual(values=ic, name="Isotype") + 
      labs(x=NULL, y=NULL) +
      theme(panel.grid.major=element_blank(), 
            panel.grid.minor=element_blank(), 
            panel.border=element_blank(),
            axis.ticks=element_blank(), 
            axis.text=element_blank(), 
            legend.position="none")
    ggsave(paste0(pa, "umap/isotype.png"))
    ggplot(subsampled, aes(umap.1, umap.2)) + 
      geom_density2d(color="#CCCCCC") +
      geom_point(size=2, aes(color=isotype)) +
      theme_bw() + 
      labs(x=NULL, y=NULL) +
      xlim(min(subsampled$umap.1)-1, max(subsampled$umap.1)+1) + 
      ylim(min(subsampled$umap.2)-1.5, max(subsampled$umap.2)+1.5) +
      scale_color_manual(values=ic, name="Isotype") + 
      labs(x=NULL, y=NULL) +
      theme(panel.grid.major=element_blank(), 
            panel.grid.minor=element_blank(), 
            panel.border=element_blank(),
            legend.text=element_text(size=15), 
            axis.ticks=element_blank(), 
            axis.text=element_blank()) + 
      guides(color=guide_legend(override.aes = list(size=7, linetype=0)))
    ggsave(paste0(pa, "umap/isotype_legend.png"))
    
    # overlay each isotype
    for (iso in isotypes) {
      ggplot(subsampled, aes(umap.1, umap.2)) + 
        geom_density2d(color="#CCCCCC") +
        geom_point(data=subsampled[isotype==iso], size=3, aes(color=iso)) +
        theme_bw() + 
        labs(x=NULL, y=NULL) +
        ggtitle(iso) +
        xlim(min(subsampled$umap.1)-1, max(subsampled$umap.1)+1) + 
        ylim(min(subsampled$umap.2)-1.5, max(subsampled$umap.2)+1.5) +
        scale_color_manual(values=isotype.colors, name="Isotype") + 
        theme(panel.grid.major=element_blank(), 
              panel.grid.minor=element_blank(), 
              panel.border=element_blank(),
              axis.ticks=element_blank(), 
              axis.text=element_blank(), 
              title=element_text(size=35),
              legend.position="none")
      ggsave(paste0(pa, "umap/isotype_overlay_", iso, ".png"))
    }
    
    for (gr in groups) {
      for (iso in isotypes) {
      ggplot(subsampled, aes(umap.1, umap.2)) + 
        geom_density2d(color="#CCCCCC") +
        geom_point(data=subsampled[group==gr & isotype==iso], size=3, aes(color=iso)) +
        theme_bw() + 
        labs(x=NULL, y=NULL) +
        xlim(min(subsampled$umap.1)-1, max(subsampled$umap.1)+1) + 
        ylim(min(subsampled$umap.2)-1.5, max(subsampled$umap.2)+1.5) +
        scale_color_manual(values=isotype.colors, name="Isotype") + 
        ggtitle(paste0(gr, "\n", iso)) +
        theme(panel.grid.major=element_blank(), 
              panel.grid.minor=element_blank(), 
              panel.border=element_blank(),
              axis.ticks=element_blank(), 
              axis.text=element_blank(), 
              title=element_text(size=35),
              legend.position="none")
      ggsave(paste0(pa, "umap/isotype_overlay_", iso, "_", gr, ".png"))
      }
    }
    
    # overlay activation molecule positive cells
    for (fac in names(yc)) {
      ggplot(subsampled, aes(umap.1, umap.2)) + 
        geom_density2d(color="#CCCCCC") +
        geom_point(data=subsampled[eval(parse(text=fac))==T], size=2, aes(color=fac)) +
        theme_bw() + 
        labs(x=NULL, y=NULL) +
        xlim(min(subsampled$umap.1)-1, max(subsampled$umap.1)+1) + 
        ylim(min(subsampled$umap.2)-1.5, max(subsampled$umap.2)+1.5) +
        ggtitle(fac) +
        #scale_color_manual(values=yc, name="") + 
        theme(panel.grid.major=element_blank(), 
              panel.grid.minor=element_blank(), 
              panel.border=element_blank(),
              axis.ticks=element_blank(), 
              axis.text=element_blank(), 
              title=element_text(size=35),
              legend.position="none")
      ggsave(paste0(pa, "umap/factors_", fac, ".png"))
    }
    
    # meta individual
    for (clust in names(mc)) {
      ggplot(subsampled, aes(umap.1, umap.2)) + 
        geom_density2d(color="#CCCCCC") +
        geom_point(data=subsampled[meta==clust], size=2, aes(color=meta)) +
        theme_bw() + 
        labs(x=NULL, y=NULL) +
        xlim(min(subsampled$umap.1)-1, max(subsampled$umap.1)+1) + 
        ylim(min(subsampled$umap.2)-1.5, max(subsampled$umap.2)+1.5) +
        labs(x=NULL, y=NULL) + 
        ggtitle(clust) +
        scale_color_manual(values=mc, name="Metacluster") + 
        theme(panel.grid.major=element_blank(), 
              panel.grid.minor=element_blank(), 
              panel.border=element_blank(),
              axis.ticks=element_blank(), 
              axis.text=element_blank(), 
              title=element_text(size=35),
              legend.position="none")
      ggsave(paste0(pa, "umap/meta_", clust, ".png"))
    }
    
    # marker
    for (marker in markers) {
      subsampled[eval(parse(text=marker))>1, eval(quote(marker)):=1]
      ggplot(subsampled, aes(umap.1, umap.2)) + 
        geom_density2d(color="#CCCCCC") +
        geom_point(aes(color=eval(parse(text=marker)), size=0.8)) +
        theme_bw() + 
        labs(x=NULL, y=NULL) +
        xlim(min(subsampled$umap.1)-1, max(subsampled$umap.1)+1) + 
        ylim(min(subsampled$umap.2)-1.5, max(subsampled$umap.2)+1.5) +
        scale_color_viridis(option="B", limits=c(-1,1)) + 
        labs(x=NULL, y=NULL) +
        ggtitle(paste(marker, "expression")) +
        theme(panel.grid.major=element_blank(), 
              panel.grid.minor=element_blank(), 
              panel.border=element_blank(),
              axis.ticks=element_blank(), 
              axis.text=element_blank(), 
              title=element_text(size=35),
              legend.position="none")
      ggsave(paste0(pa, "umap/", marker, ".png"))
    }
    
    # overlay age group
    for (fa in groups) {
      ggplot(subsampled, aes(umap.1, umap.2)) + 
        geom_density2d(color="#CCCCCC") +
        geom_point(data=subsampled[group==fa], size=3, aes(color=fa)) +
        theme_bw() + 
        labs(x=NULL, y=NULL) +
        xlim(min(subsampled$umap.1)-1, max(subsampled$umap.1)+1) + 
        ylim(min(subsampled$umap.2)-1.5, max(subsampled$umap.2)+1.5) +
        scale_color_manual(values=group.colors, name="Group") + 
        ggtitle(fa) +
        theme(panel.grid.major=element_blank(), 
              panel.grid.minor=element_blank(), 
              panel.border=element_blank(),
              axis.ticks=element_blank(), 
              axis.text=element_blank(), 
              title=element_text(size=35),
              legend.position="none")
      ggsave(paste0(pa, "umap/group_overlay_", fa, ".png"))
    }
    
    # overlay age group
    # day 7 only
    for (fa in groups) {
      ggplot(subsampled, aes(umap.1, umap.2)) + 
        geom_density2d(color="#CCCCCC") +
        geom_point(data=subsampled[group==fa & timepoint=="day7"], size=3, aes(color=fa)) +
        theme_bw() + 
        labs(x=NULL, y=NULL) +
        xlim(min(subsampled$umap.1)-1, max(subsampled$umap.1)+1) + 
        ylim(min(subsampled$umap.2)-1.5, max(subsampled$umap.2)+1.5) +
        scale_color_manual(values=group.colors, name="Group") + 
        ggtitle(paste0(fa, "\nDay 7 Total B cells")) +        
        theme(panel.grid.major=element_blank(), 
              panel.grid.minor=element_blank(), 
              panel.border=element_blank(),
              axis.ticks=element_blank(), 
              axis.text=element_blank(), 
              title=element_text(size=35),
              legend.position="none")
      ggsave(paste0(pa, "umap/day7_group_overlay_", fa, ".png"))
    }
    
    subsampled$density <- get_density(subsampled$umap.1, subsampled$umap.2, n=500)
    for (fa in groups) {
      ggplot(subsampled, aes(umap.1, umap.2)) + 
        geom_density2d(color="#CCCCCC") +
        geom_point(data=subsampled[group==fa & timepoint=="day7"], size=3, aes(color=density)) +
        theme_bw() + 
        labs(x=NULL, y=NULL) +
        xlim(min(subsampled$umap.1)-1, max(subsampled$umap.1)+1) + 
        ylim(min(subsampled$umap.2)-1.5, max(subsampled$umap.2)+1.5) +
        scale_color_viridis(option='magma', name="Density", direction=1) + 
        ggtitle(paste0(fa, "\nDay 7 Total B cells")) +
        theme(panel.grid.major=element_blank(), 
              panel.grid.minor=element_blank(), 
              panel.border=element_blank(),
              axis.ticks=element_blank(), 
              axis.text=element_blank(), 
              title=element_text(size=35),
              legend.position="none")
      ggsave(paste0(pa, "umap/day7_density_group_overlay_", fa, ".png"))
    }
    
  }
