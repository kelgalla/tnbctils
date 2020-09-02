library(ggplot2)
gn <- c(rep("gn1", 5),  rep("gn2", 5),  rep("gn3", 5),  rep("gn4", 5),  rep("gn5", 5), 
        rep("gn6", 5),  rep("gn7", 5),  rep("gn8", 5),  rep("gn9", 5),  rep("gn10", 5),
        rep("gn11", 5), rep("gn12", 5), rep("gn13", 5), rep("gn14", 5), rep("gn15", 5),
        rep("gn16", 5), rep("gn17", 5), rep("gn18", 5), rep("gn19", 5))

smpl <- rep(c("smpl1", "smpl2", "smpl3", "smpl4", "smpl5"), 19)

mut1 <- c(0, 2, 3, 3, 0)
mut2 <- c(0, 2, 2, 0, 3)
mut3 <- c(2, 0, 2, 0, 3)
mut4 <- c(0, 0, 2, 2, 3)
mut5 <- c(0, 2, 3, 3, 2)

mut <- c(mut1, mut2, mut3, mut4, mut5, mut2, mut3, mut1, mut4, mut2,
         mut5, mut1, mut4, mut5, mut5, mut1, mut3, mut1, mut4)

mut <- factor(mut, levels=c(0,2,3))

tst <- data.frame(gn=gn, smpl=smpl, mut=mut)
tst <- tst[1:25,]

n1 <- length(unique(tst$gn))
n2 <- length(unique(tst$smpl))

png(paste0("F:\\TNBC TILS\\tst.png"))

ggplot(tst) + 
  geom_tile(aes(x=smpl, y=gn, fill=mut)) +
  guides(fill=F) +
  ylab("Genes") +
  xlab("Sample") +
  scale_x_discrete(expand = expand_scale(add =0)) +
  scale_y_discrete(expand = expand_scale(add=0)) +
  theme(panel.border=element_rect(colour="#000000", fill=NA)) +
  scale_fill_manual(values=c("#FFFFFF", "#00BFC4", "#F8766D")) +
  theme(panel.grid.major.y=element_blank()) +
  theme(panel.grid.major.x=element_blank()) +
  theme(panel.grid.minor.y=element_blank()) +
  theme(panel.grid.minor.x=element_blank()) +
  theme(panel.background=element_blank()) +
  theme(axis.text.y = element_text(size=20, colour="black", face="bold")) +
  theme(axis.title = element_text(size=22, colour="black", face="bold")) +
  theme(plot.title = element_text(size=24, colour="black", face="bold")) +
  theme(axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.x=element_blank())
  #geom_line(data = data.frame(x = c(0, n2) + 0.5, y = rep(2:n1, each = 2) - 0.5),
           # aes(x = x, y = y, group = y))

dev.off()

