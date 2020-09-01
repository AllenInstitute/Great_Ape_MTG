library(pheatmap)
library(ggplot2)

calc_auc <- function(x, y) {
  f <- approxfun(data.frame(x, y))
  integrate(f, min(x), max(x))$value
}


# Load overlaps
f.path <- "//allen/programs/celltypes/workgroups/rnaseqanalysis/Nik/Projects_for_other_people/Trygve_Bakken/cluster_overlap_metric_test_data/"
fn.l <- list.files(f.path, pattern = "*.csv")

m.l <- list()
for (fn1 in fn.l) {
  res1 <- sub(".csv", "", strsplit(fn1, "_")[[1]][3], fixed = TRUE)
  m1 <- read.csv(file = paste0(f.path, fn1), row.names = 1)
  m.l[[res1]] <- as.matrix(m1)
}
res.order <- order(as.numeric(sub("res", "", names(m.l))))
m.l <- m.l[res.order]

m1 <- m.l[[1]]
pheatmap(m1, cluster_rows = FALSE, cluster_cols = FALSE)



min.overlap <- 0.1  # Min cluster overlap to be considered a match
th.step <- 0.1
prop.th <- seq(min.overlap, 1, by = th.step)

auc.df <- data.frame()
prec.recall.df <- data.frame()
homol <- NULL
m1.bin.l <- list()
for (res1 in names(m.l)) {
  m1 <- m.l[[res1]]
  ms <- list(m1, t(m1))
  
  m1.bin.res <- m1 * 0
  auc.m <- NULL
  prec.max <- NULL
  recall.max <- NULL
  for (k in 1:2) {
    m1 <- ms[[k]]
    prec <- NULL
    recall <- NULL
    
    for (th1 in prop.th) {
      m1.bin <- ifelse(m1 > th1, th.step, 0)
      
      cl.match.cnt <- apply(m1, 1, function(x) sum(x > th1))
      recall1 <- sum(cl.match.cnt > 0) / nrow(m1)
      if (recall1 > 0) {
        recall <- c(recall, recall1)
        res <- mean(cl.match.cnt[cl.match.cnt > 0])
        prec1 <- (ncol(m1) / res - 1) / (ncol(m1) - 1)
        prec <- c(prec, prec1)
        prec.recall.df <- rbind(prec.recall.df, 
                                data.frame(cl_resolution = res1, overlap_th = th1, 
                                           dataset = as.factor(k), 
                                           precision = prec1, recall = recall1))
        
        # Define homologies
        cl.match.col <- apply(m1, 1, function(x) names(which(x > th1)))
        cl.match.row <- apply(m1, 2, function(x) names(which(x > th1)))
        for (i in 1:length(cl.match.col)) {
          col1 <- cl.match.col[[i]]
          for (j in 1:10) {
            row1 <- unique(unlist(cl.match.row[col1]))
            col1 <- unique(unlist(cl.match.col[row1]))
          }
          m1.bin[row1, col1] <- th.step
        }
        if (k == 1) m1.bin.res <- m1.bin.res + m1.bin
      }
    }
    m1.bin.l[[res1]] <- m1.bin.res
    
    prec.max <- c(prec.max, max(prec))
    recall.max <- c(recall.max, max(recall))
    auc1 <- calc_auc(c(0, recall), c(0, prec))
    auc.m <- c(auc.m, auc1)
  }
  # plot(recall, prec, xlim = c(0, 1), ylim = c(0, 1), main = res1)
  auc.df <- rbind(auc.df, data.frame(res = res1, auc = mean(auc.m), 
                                     max_prec = min(prec.max),
                                     max_recall = min(recall.max)))
}

plot(auc.df$max_recall, auc.df$auc, ylim = c(0, 1))
text(auc.df$max_recall, auc.df$auc, labels = auc.df$res, pos = 3)

g.th <- ggplot(prec.recall.df, aes(x = recall, y = precision, 
                                   color = overlap_th, shape = dataset)) +
  facet_wrap(~ cl_resolution) +
  geom_point() +
  theme_bw()
plot(g.th)





# Identify homologies
m2 <- m.l[[4]]
m2.bin <- ifelse(m1.bin.l[[4]] > 0.4, 1, 0)
pheatmap(m2.bin, cluster_rows = FALSE, cluster_cols = FALSE)

library(biclust)
m2.bc <- biclust(m2.bin, method = BCBimax(), minr=1, minc=1)
m2.bc.memb <- biclusternumber(m2.bc)
heatmapBC(m2, m2.bc)

# m1.bc.l[[res1]][[as.character(th1)]] <- m2.bc.memb






#### Define homologies ####

m1 <- m.l[[1]]

th1 <- 0.25
m1.bin <- ifelse(m1 > th1, 1, 0)
pheatmap(m1, cluster_rows = FALSE, cluster_cols = FALSE)
pheatmap(m1.bin, cluster_rows = FALSE, cluster_cols = FALSE)

cl.match.col <- apply(m1, 1, function(x) names(which(x > th1)))
cl.match.row <- apply(m1, 2, function(x) names(which(x > th1)))
for (i in names(cl.match.col)) {
  col1 <- cl.match.col[[i]]
  for (i in 1:10) {
    row1 <- unique(unlist(cl.match.row[col1]))
    col1 <- unique(unlist(cl.match.col[row1]))
  }
  row.col <- update_cl(cl.match.col, cl.match.row, col1)
  m1.bin[row.col[[1]], row.col[[2]]] <- 1
}

m1.bc <- biclust(m1.bin, method = BCBimax(), minr=1, minc=1)
# bc2 <- biclust(m1, method = BCCC())

# bicluster(m1.bin, m1.bc)
biclusternumber(m1.bc)
