library(ggsci)
library(ggplot2)
library(grid)
library(gridExtra)
library(rjson)
library(Cairo)
library(cowplot)

names_of_sets = c(" full NCBI anno.", " detected 6frame", " detected NCBI anno.")
#parameters <- fromJSON(file = "./parameters.json")
parameters <- fromJSON(file = "~/Documents/paper/PROTMAP_clean/scripts/parameters.json")

sf1 <- paste(parameters$data_dir, "/accumulated_data/all_CDS", sep = "")
set1 <- array(read.table(sf1)[,1])
Col <- pal_npg("nrc", alpha=0.8)(10)
Col <- c(Col[1], Col[4], Col[3], Col[9], Col[6], Col[7])

names <- c("NCBI anno", "6frame", "NCBI","both", "only 6frame","only NCBI")

k <- 1
sf2 <- paste(parameters$data_dir, "/accumulated_data/k_", k, "_6frame", sep = "")
sf3 <- paste(parameters$data_dir, "/accumulated_data/k_", k, "_proteom", sep = "")
set2 <- array(read.table(sf2)[,1])
set3 <- array(read.table(sf3)[,1])
df <- data.frame(sets=names,
                 PSMs=c(length(set1), length(set3), length(set2),
                        length(intersect(set2, set3)), length(setdiff(set3, set2)),
                        length(setdiff(set2, set3))))
df$sets <- factor(df$sets,levels=names)
p <- ggplot(data = df, aes(x=sets, y=PSMs, fill=sets)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = Col) +
  theme(legend.position="bottom", legend.direction = "horizontal",
        legend.title = element_text(size = 0),
        legend.text = element_text(size = 45),
        legend.spacing.x = unit(1, "lines"))
legend <- cowplot::get_legend(p)
pdf("~/Documents/paper/PROTMAP_clean/figs/data_base_compare/barplot_legend.pdf", width = 11, height = 3) 
#CairoPDF("../figs/data_base_compare/barplot_legend.pdf", width = 7, height = 7)
grid.newpage()
grid.draw(legend)
dev.off()
p <- ggplot(data = df, aes(x=sets, y=PSMs, fill=sets)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = Col) +
  coord_flip() +
  annotate("text", x=1, y=3100, label="4000", size=6) +
  theme(legend.position="none",
        axis.title.y = element_blank(),
        axis.text = element_text(size = 25),
        axis.ticks.y = element_blank(),
        axis.title = element_text(size = 35))

pdf(paste("~/Documents/paper/PROTMAP_clean/figs/data_base_compare/barplot_", k, ".pdf", sep = ""), width = 9, height = 7) 
#CairoPDF(paste("../figs/data_base_compare/venn_", k, ".pdf", sep = ""), width = 7, height = 7)
plot(p)
dev.off()

for(k in c(6,10)){
  sf2 <- paste(parameters$data_dir, "/accumulated_data/k_", k, "_6frame", sep = "")
  sf3 <- paste(parameters$data_dir, "/accumulated_data/k_", k, "_proteom", sep = "")
  set2 <- array(read.table(sf2)[,1])
  set3 <- array(read.table(sf3)[,1])
  df <- data.frame(sets=names,
                   PSMs=c(length(set1), length(set3), length(set2),
                          length(intersect(set2, set3)), length(setdiff(set3, set2)),
                          length(setdiff(set2, set3))))
  df$sets <- factor(df$sets,levels=names)
  p <- ggplot(data = df, aes(x=sets, y=PSMs, fill=sets)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = Col) +
    coord_flip() +
    theme(legend.position="none",
          axis.title.y = element_blank(),
          axis.text = element_text(size = 25),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title = element_text(size = 35))
  pdf(paste("~/Documents/paper/PROTMAP_clean/figs/data_base_compare/barplot_", k, ".pdf", sep = ""), width = 7, height = 7) 
  #CairoPDF(paste("../figs/data_base_compare/venn_", k, ".pdf", sep = ""), width = 7, height = 7)
  plot(p)
  dev.off()
}

for(k in c(1,6,10)){
  sf2 <- paste(parameters$data_dir, "/accumulated_data/k_", k, "_6frame", sep = "")
  sf3 <- paste(parameters$data_dir, "/accumulated_data/k_", k, "_proteom", sep = "")
  set2 <- array(read.table(sf2)[,1])
  set3 <- array(read.table(sf3)[,1])
  df <- data.frame(sets=names,
                   PSMs=c(length(set1), length(set3), length(set2),
                          length(intersect(set2, set3)), length(setdiff(set3, set2)),
                          length(setdiff(set2, set3))))
  df$sets <- factor(df$sets,levels=names)
  p <- ggplot(data = df, aes(x=sets, y=PSMs, fill=sets)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = Col) +
    coord_flip() +
    theme(legend.position="none",
          axis.title.y = element_blank(),
          axis.text = element_text(size = 25),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title = element_text(size = 35))
  pdf(paste("~/Documents/paper/PROTMAP_clean/figs/data_base_compare/barplot_", k, ".pdf", sep = ""), width = 7, height = 7) 
  #CairoPDF(paste("../figs/data_base_compare/venn_", k, ".pdf", sep = ""), width = 7, height = 7)
  plot(p)
  dev.off()
}

PSM <- c(rep("1" , 5) , rep("6" , 5), rep("10" , 5))
sub_set_names <- c("6frame not in NCBI", "only 6frame", "both", "only NCBI", "NCBI not called")
sub_sets <- rep(sub_set_names , 3)
full_NCBI <- paste(parameters$data_dir, "/accumulated_data/all_CDS", sep = "")
full_NCBI <- array(read.table(full_NCBI)[,1])
values <- numeric()
for(k in c(1,6,10)){
  frame <- paste(parameters$data_dir, "/accumulated_data/k_", k, "_6frame", sep = "")
  NCBI <- paste(parameters$data_dir, "/accumulated_data/k_", k, "_proteom", sep = "")
  frame <- array(read.table(frame)[,1])
  NCBI <- array(read.table(NCBI)[,1])
  
  frame_not_in_anno <- length(setdiff(frame, full_NCBI))
  only_frame_in_anno <- length(setdiff(intersect(frame, full_NCBI), NCBI))
  both <- length(intersect(frame, NCBI))
  only_NCBI <- length(setdiff(NCBI, frame))
  NCBI_not_called <- length(setdiff(setdiff(full_NCBI, frame), NCBI))
  values <- c(values, frame_not_in_anno, only_frame_in_anno, both, only_NCBI, NCBI_not_called)
}

Col <- pal_npg("nrc", alpha=0.7)(10)
Col <- c(Col[8], Col[1], Col[9], Col[3], Col[2])

data <- data.frame(PSM,sub_sets,values)
data$sub_sets <- factor(data$sub_sets, levels=sub_set_names)
data$PSM <- factor(data$PSM, levels=c("1", "6", "10"))
p <- ggplot(data, aes(fill=sub_sets, y=values, x=PSM)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = Col) + 
  geom_segment(x=0.5, xend=3.5, y=length(full_NCBI), yend=length(full_NCBI),
               linetype="dashed", color = "black", size=0.5) +
  labs(y = "Proteins", fill = "") + 
  guides(fill=guide_legend(nrow=2,byrow=FALSE)) +
  theme(legend.position="none",
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 35),
        legend.text = element_text(size = 25))

pdf(paste("~/Documents/paper/PROTMAP_clean/figs/data_base_compare/barplot_all.pdf", sep = ""), width = 10, height = 7) 
plot(p)
dev.off()



#CairoPDF("../figs/data_base_compare/barplot_legend.pdf", width = 7, height = 7)
PSM <- c(rep("1" , 4))
sub_set_names <- c("only 6frame", "both", "only NCBI", "NCBI not called")
values <- c(1,2,3,4)
data <- data.frame(PSM,sub_set_names,values)
data$sub_set_names <- factor(data$sub_set_names, levels=sub_set_names)
Col <- pal_npg("nrc", alpha=0.7)(10)
Col <- c(Col[8], Col[9], Col[3], Col[2])
p <- ggplot(data, aes(fill=sub_set_names, y=values, x=PSM)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = Col) + 
  guides(fill=guide_legend(nrow=2,byrow=FALSE)) +
  theme(legend.position="bottom",
        legend.title = element_text(size = 0),
        legend.spacing.x = unit(0.5, 'cm'))
pdf("~/Documents/paper/PROTMAP_clean/figs/data_base_compare/barplot_legend.pdf", width = 3, height = 1) 
legend <- cowplot::get_legend(p)
grid.newpage()
grid.draw(legend)
dev.off()
plot(p)

PSM <- c(rep("1" , 4) , rep("6" , 4), rep("10" , 4))
sub_set_names <- c("only 6frame", "both", "only NCBI", "NCBI not called")
sub_sets <- rep(sub_set_names , 3)
full_NCBI <- paste(parameters$data_dir, "/accumulated_data/all_CDS", sep = "")
full_NCBI <- array(read.table(full_NCBI)[,1])
values <- numeric()
for(k in c(1,6,10)){
  frame <- paste(parameters$data_dir, "/accumulated_data/k_", k, "_6frame", sep = "")
  NCBI <- paste(parameters$data_dir, "/accumulated_data/k_", k, "_proteom", sep = "")
  frame <- array(read.table(frame)[,1])
  NCBI <- array(read.table(NCBI)[,1])
  
  only_frame <-length(setdiff(frame, NCBI))
  both <- length(intersect(frame, NCBI))
  only_NCBI <- length(setdiff(NCBI, frame))
  NCBI_not_called <- length(setdiff(setdiff(full_NCBI, frame), NCBI))
  values <- c(values, only_frame, both, only_NCBI, NCBI_not_called)
}

Col <- pal_npg("nrc", alpha=0.7)(10)
Col <- c(Col[8], Col[9], Col[3], Col[2])

data <- data.frame(PSM,sub_sets,values)
data$sub_sets <- factor(data$sub_sets, levels=sub_set_names)
data$PSM <- factor(data$PSM, levels=c("1", "6", "10"))
p <- ggplot(data, aes(fill=sub_sets, y=values, x=PSM)) + 
  geom_bar(position="stack", stat="identity",key_glyph = "polygon3") +
  scale_fill_manual(values = Col) + 
  geom_segment(x=0.5, xend=3.5, y=length(full_NCBI), yend=length(full_NCBI),
               linetype="dashed", color = "black", size=0.5) +
  labs(y = "Proteins", fill = "") + 
  guides(fill=guide_legend(nrow=2, byrow=FALSE)) +
  theme(legend.position="right",
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 35),
        legend.spacing = unit(1.0, "cm"),
        legend.key.size = unit(4, "cm"),
        legend.margin = margin(0,0,0,200),
        legend.text = element_text(size = 35),
        legend.background=element_blank(),
        legend.key = element_rect(colour = "transparent", fill = "white"))+
pdf(paste("~/Documents/paper/PROTMAP_clean/figs/data_base_compare/barplot_all.pdf", sep = ""), width = 22, height = 7) 
print(p)
dev.off()
draw_key_polygon3 <- function(data, params, size) {
  lwd <- min(data$size, min(size) / 4)
  
  grid::rectGrob(
    width = grid::unit(0.5, "npc"),
    height = grid::unit(0.6, "npc"),
    gp = grid::gpar(
      col = data$colour,
      fill = alpha(data$fill, data$alpha),
      lty = data$linetype,
      lwd = lwd * .pt,
      linejoin = "mitre"
    ))
}
GeomBar$draw_key = draw_key_polygon3

### this step is not needed anymore per tjebo's comment below
### see also: https://ggplot2.tidyverse.org/reference/draw_key.html
# register new key drawing function, 
# the effect is global & persistent throughout the R session
