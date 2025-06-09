rm(list = ls())

library(VennDiagram)
library(grid)

total_genes <- 1180
te_early <- 162
totalRNA_late <- 294
overlap <- 82

p_value <- phyper(overlap - 1, te_early, total_genes - te_early, totalRNA_late, lower.tail = FALSE)
enrichment_ratio <- (overlap / totalRNA_late) / (te_early / total_genes)

stats_text <- paste(
  sprintf("Total 1180 genes"),
  sprintf("Hypergeometric test p-value = %.2e", p_value),
  sep = "\n"
)

venn.plot <- draw.pairwise.venn(
  area1 = totalRNA_late,
  area2 = te_early,
  cross.area = overlap,
  category = c("totalRNA late upregulated", "TE early upregulated"),
  fill = c("#e06666", "#6287e0"),
  lty = "blank",
  cex = 2.0,
  cat.cex = 1.8,
  cat.pos = c(-20, 20), 
  cat.dist = 0.05,
  scaled = TRUE
)

grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 1, heights = unit(c(1, 8), "null"))))

grid.text(stats_text, vp = viewport(layout.pos.row = 1, layout.pos.col = 1), gp = gpar(fontsize = 25))

grid.draw(gTree(children = venn.plot, vp = viewport(layout.pos.row = 2, layout.pos.col = 1)))

png("venn_with_stats.png", width = 600, height = 700)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 1, heights = unit(c(1, 8), "null"))))
grid.text(stats_text, vp = viewport(layout.pos.row = 1, layout.pos.col = 1), gp = gpar(fontsize = 12))
grid.draw(gTree(children = venn.plot, vp = viewport(layout.pos.row = 2, layout.pos.col = 1)))
dev.off()