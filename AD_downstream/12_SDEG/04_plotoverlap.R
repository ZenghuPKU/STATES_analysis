rm(list = ls())

library(VennDiagram)
library(grid)

total_genes <- 2837
SDEG <- 66
PIG <- 57
overlap <- 27

p_value <- phyper(overlap - 1, SDEG, total_genes - SDEG, PIG, lower.tail = FALSE)
enrichment_ratio <- (overlap / PIG) / (SDEG / total_genes)

stats_text <- paste(
  sprintf("Total %d genes", total_genes),
  sprintf("Hypergeometric test p-value = %.2e", p_value),
  sep = "\n"
)

venn.plot <- draw.pairwise.venn(
  area1 = PIG,
  area2 = SDEG,
  cross.area = overlap,
  category = c("PIG", "SDEG"),
  fill = c("#e06666", "#6287e0"),
  lty = "blank",
  cex = 2.0,
  cat.cex = 1.8,
  cat.pos = c(-20, 20), 
  cat.dist = 0.05,
  scaled = TRUE
)

pdf("venn_with_stats.pdf", width = 6, height = 7)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 1, heights = unit(c(1, 8), "null"))))
grid.text(stats_text, vp = viewport(layout.pos.row = 1, layout.pos.col = 1), gp = gpar(fontsize = 18))
grid.draw(gTree(children = venn.plot, vp = viewport(layout.pos.row = 2, layout.pos.col = 1)))
dev.off()