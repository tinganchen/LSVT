library(readxl)
# Data input with the 128TH var.(data length) skipped
data <- read_excel("data\\LSVT_voice_rehabilitation.xlsx",
                   sheet = 1, col_names = T)[ , -128]
d.class <- read_excel("data\\LSVT_voice_rehabilitation.xlsx",
                      sheet = 2, col_names = T)
d.subject <- read_excel("data\\LSVT_voice_rehabilitation.xlsx",
                        sheet = 3, col_names = T)
colnames(d.class) <- "class"
colnames(d.subject) <- c("id", "age", "sex")
class <- d.class$class
dim(data)
dim(d.class)
dim(d.subject)


#============================================================

# EDA
## Subject
# gender
pie.lab <- paste(c("Male", "Female"), "\n", table(d.subject$sex) / 9)
pie(table(d.subject$sex) / 9, labels = pie.lab,
    col = c("blue", "red"), clockwise = T,
    main = "PD subjects", cex = 1.1, density  = 30)

# age, gender
d.age.sex <- data.frame(table(d.subject[ , -1]) / 9)
barplot(t(cbind(d.age.sex[1:10 , 3], d.age.sex[11:20 , 3])),
        beside = T, names.arg = d.age.sex[1:10, 1],
        main = "PD subjects", ylim = c(0, 3),
        xlab = "Age", ylab = "Number of Male/Female",
        col = c("cornflowerblue", "pink"),
        legend.text = c("Male", "Female"),
        args.legend = "topright")

## class(whether the treatment effect is good(acceptable))
prop.phon.judge <- round(prop.table(table(d.class)) * 100)

pie.lab2 <- paste(c("Acceptable", "Unacceptable"), "\n", 
                  table(d.class), "(", prop.phon.judge, "%)")
pie(table(d.class), labels = pie.lab2,
    col = c("blue", "red"), clockwise = T,
    main = "Phonation", cex = 1.1, density  = 30)


## Relationship between algorithm & class(treatment effect)
# correlations 
corr.algo.cl <- apply(data, 2, function(x) cor(x, class))
top10.corr.var <- head(names(sort(abs(corr.algo.cl),
                                  decreasing = T)), 10)
top10.high.corr <- corr.algo.cl[top10.corr.var]
last10.corr.var <- head(names(sort(abs(corr.algo.cl))), 10)
last10.corr <- corr.algo.cl[last10.corr.var]

## Density plot 
library(sm)
# accept./ unaccept.
class.f <- factor(class, levels = c(1, 2),
                  labels = c("Acceptable", "Unacceptable"))
denplot.fun <- function(algor, cor) {
  var <- as.numeric(as.matrix(data[ , algor]))
  sm.density.compare(var, class.f, xlab = "", 
                     col = c(4, 2), lwd = 2)
  title(paste("\n\n", algor))
  legend("center", bty = "n", cex = 1.5,
         paste("Correlation", round(cor[algor], 2)))
}
par(mfrow = c(3, 4))
par(mai = c(0.4, 0.6, 0.6, 0.2))
# Top10
top10.denplot <- apply(as.matrix(top10.corr.var), 1, 
                       denplot.fun, top10.high.corr)
plot(0, xaxt = "n", yaxt = "n", type = "n",
     xlab = "", ylab = "", main = "\n\nJudgment")
legend("center", levels(class.f), bty = "n", cex = 1.2,
       col = c(4, 2), lty = 1:2, lwd = 2)
title("\nTop10 correlated algorithm Density Estimates", 
      outer = T, cex.main = 2)
# Last10
par(mfrow = c(3, 4))
par(mai = c(0.4, 0.6, 0.6, 0.2))
last10.denplot <- apply(as.matrix(last10.corr.var), 1, 
                        denplot.fun, last10.corr)
plot(0, xaxt = "n", yaxt = "n", type = "n",
     xlab = "", ylab = "", main = "Judgment")
legend("center", levels(class.f), bty = "n", cex = 1.2,
       col = c(4, 2), lty = 1:2, lwd = 2)
title("\n10 Least correlated algorithm Density Estimates", 
      outer = T, cex.main = 2)

# Gender/object
par(mfrow = c(3, 4))
par(mai = c(0.6, 0.6, 0.5, 0.2))
apply(as.matrix(top10.corr.var), 1, function(x) {
  var <- as.numeric(as.matrix(data[ , x]))
  sm.density.compare(var, d.subject$sex, xlab = x, 
                     col = c(4, 2))
})
plot(0, xaxt = "n", yaxt = "n", type = "n",
     xlab = "", ylab = "")
legend("center", c("Male", "Female"), bty = "n",
       col = c(4, 2), lty = 1:2, title = "Gender")
title("\n\nDensity Estimates", outer = T)

# Age/object

#=================================================

## Feature & Dimension Reduction method Selection
# (criterion: Max. classification(LDA) accuracy)

# 1. Data updated by 
#    Oredered significant variables (algorithm)
bw.ratio <- function(x, y){
  tg <- table(y)
  gm <- tapply(x, y, mean)
  bw <- (sum((gm-mean(x))^2))/(sum((x - rep(gm, tg))^2))
}
bw.values <- apply(data, 2, bw.ratio, class)
order.bw <- order(bw.values, decreasing = T)
d.algor.ord <- as.matrix(data)[ , order.bw]

# 2. LDA, leave-one-out, function
library(MASS)
lda.accur.fun <- function(dim, data) {
  lda.data <- lda(x = data[ , 1:dim], 
                  grouping = class.f, CV = T)
  # In case of NaN
  keep <- unique(unlist(apply(lda.data$posterior, 2, function(x)
    which(is.nan(x) != 1))))
  confusion.mat <- table(class[keep], 
                         lda.data$class[keep])
  accuracy <- sum(diag(prop.table(confusion.mat))) 
  accuracy
}

# 2. Dimension Reduction method Selection (matplot.)
algor.DR.method.select.fun <- function(algor.num) {
  LSVT <- d.algor.ord[ , 1:algor.num]
  # a). PC scores
  library(FactoMineR)
  library(ggplot2)
  library(factoextra)
  lsvt.pca <- PCA(LSVT, scale.unit = T, graph = F, 
                  ncp = dim(LSVT)[1] - 1)
  eig <- get_eigenvalue(lsvt.pca)
  var <- get_pca_var(lsvt.pca)
  ind <- get_pca_ind(lsvt.pca)
  lsvt.pc.score <- ind$coord[ , 1:10]
  lda.pc.score.accur <- apply(as.matrix(2:10), 1, 
                              lda.accur.fun, lsvt.pc.score)
  # b). MDS projected value
  mds.subj <- cmdscale(dist(scale(LSVT)), k = 10)
  lda.mda.accur <- apply(as.matrix(2:10), 1, 
                         lda.accur.fun, mds.subj)
  # c). ISOMAP projected value
  library(vegan)
  iso.score.k <- apply(as.matrix(2:20), 1, function(neighbor) {
    lsvt.isomap <- isomap(dist(scale(LSVT)), 
                          ndim = 10, k = neighbor)$points
    dim.accur <- apply(as.matrix(2:10), 1, 
                       lda.accur.fun, lsvt.isomap)
    dim.accur
  })
  lda.iso.accur <- apply(iso.score.k, 1, max)
  
  lda.accur <- rbind(rep(NA, 3),
                     data.frame(PCA = lda.pc.score.accur,
                                MDS = lda.mda.accur,
                                ISOMAP = lda.iso.accur))
  matplot(1 - lda.accur, pch = c(1, 3, 2), type = "b",
          col = 2:4, lwd = 2, cex = c(2, 1.5, 1.2),
          main = paste("\n", algor.num, " algorithms"),
          xlab = "Dimensions", ylab = "Error rate")
}
# Before feature selection
algor.DR.method.select.fun(dim(d.algor.ord)[2])
legend("topright", c("PCA", "MDS", "ISOMAP"),
       pch = c(1, 3, 2), lwd = 2, col = 2:4,
       lty = 2:4, pt.cex = c(2, 1.5, 1.2))
# After feature selection (ISOMAP is always the best)
par(mfrow = c(3, 3))
par(mai = c(0.6, 0.6, 0.5, 0.2))
algor.DR.method.select <- apply(as.matrix(c(20, 40, 60, 
                                            80, 100, 150, 
                                            200, 250, 309)),
                                1, algor.DR.method.select.fun)
legend("top", c("PCA", "MDS", "ISOMAP"),
       pch = c(1, 3, 2), lwd = 2, col = 2:4, horiz = T,
       lty = 1:3, pt.cex = c(1.5, 1.2, 1.2), cex = 0.8,
       bg = cm.colors(1, alpha = 0), bty = "n")
title("\nError Rate by LDA upon Variables Selection",
      outer = T)

# 3. Feature selection
# Compare isomap error rate in diff. number of var. selected
algor.isomap.select.fun <- function(algor.num) {
  LSVT <- d.algor.ord[ , 1:algor.num]
  # c). ISOMAP projected value
  iso.score.k <- apply(as.matrix(2:20), 1, function(neighbor) {
    lsvt.isomap <- isomap(dist(scale(LSVT)), 
                          ndim = 10, k = neighbor)$points
    dim.accur <- apply(as.matrix(2:10), 1, 
                       lda.accur.fun, lsvt.isomap)
    dim.accur
  })
  lda.iso.accur <- apply(iso.score.k, 1, max)
  lda.iso.accur
}
algor.isomap.select <- rbind(rep(NA, 9),
                             apply(as.matrix(c(20, 40, 60, 
                                               80, 100, 150, 
                                               200, 250, 309)),
                                   1, algor.isomap.select.fun))
library(fields)
matplot(1-algor.isomap.select, pch = 16,
        type = "l", lwd = 4, lty = 1,
        xlim = c(0.5, 10.5), ylim = c(0.06, 0.28),
        col = tim.colors(9), xlab = "Dimensions",
        ylab = "Error rate", 
        main = "ISOMAP Classification Error Rate
        (Leave-one-out LDA)",
        xaxt = "n")
axis(1, at = 2:10, 2:10)
legend("topright", legend = rev(c(20, 40, 60, 80, 100, 150, 
                                  200, 250, 309)),
       lwd = 3, col = rev(tim.colors(9)), ncol = 3,
       title = "Feature selection", cex = 0.8)
points(rep(7, 2), 
       rep(1-algor.isomap.select[8, 2], 2),
       pch = c(16, 1), col = 2, cex = c(1, 3)) # (Continue...)
# Thus, the first 40 sig. var.(features) are selected
# to classify trt. consequences (LDA) 
# after Dimension Reduction (ISOMAP dim.7).

# 4. Data that contains only the 40 most sig. var.
data.sig.var <- d.algor.ord[ , 1:40]

# 5. ISOMAP dim.7, LDA accuracy
##   Isomap dim.1-7
# a). Find the most proper k(# of neighbors)
accur.given.k <- apply(as.matrix(2:20), 1, function(x) {
  iso.pts <- isomap(dist(scale(data.sig.var)), 
                    ndim = 7, k = x)$points
  accuracy <- lda.accur.fun(dim = 7, data = iso.pts)
  accuracy 
})
best.k <- (2:20)[which.max(accur.given.k)]

# b). Isomap
lsvt.isomap <- isomap(dist(scale(data.sig.var)), 
                      ndim = 7, k = best.k)
# c). classification by LDA
lda.isomap <- lda(x = lsvt.isomap$points, 
                  grouping = class.f, CV = T)
iso.best.accur <- sum(diag(prop.table(
  table(class, lda.isomap$class)))) # 91.27%

# (continue...accuracy written in)
text(6, 1-algor.isomap.select[8, 2]-0.01,
     paste("Accuracy\n= ", round(100*iso.best.accur), "%"),
     col = 2, cex = 1.2)
# Mean error rate in Different numbers of var. selected
mean.err.rate <- 1-apply(algor.isomap.select[-1, ], 2, mean)
par(new = T)
plot(rep(1, 9), mean.err.rate, 
     col = tim.colors(9), pch = 16, cex = 1.5,
     xlab = "", ylab = "", xaxt = "n", yaxt = "n",
     bty = "n", xlim = c(0.5, 10.5), ylim = c(0.06, 0.28))
axis(1, at = 1, "Mean\nerror rate", tick = "n")


# 6. Find the cause of accept./unaccept. of phonation
#  class(LDA posterior) => dim.val.(isomap$points)
#  => algorithm [measured by correlation coefficients]

# a). class(LDA posterior) => dim.val.(isomap$points)
#    Effects of isomap$points (dim.) 
#    on the class (LSVT treatment effect)
# Measured by the corr. between isomap$points &
# posterior prob.(the prob. that the phonation is
# accept./unaccept.)
cor.post.dim <- apply(lsvt.isomap$points, 2, function(x)
  cor(x, lda.isomap$posterior[ , 2]))
ord.dim.eff <- order(abs(cor.post.dim), decreasing = T)
reord.cor <- cor.iso.dim.post[ord.dim.eff]

library(fields)

col.red.blue <- two.colors(start = "blue", 
                           middle = "gray", 
                           end = "red")
trans <- function(x) {round((x+1)*127.5)+1}

plot(1:7, reord.cor, cex = 20 * (abs(reord.cor)),
     xaxt = "n", xlim = c(0.5, 7.5), ylim = c(-1, 1),
     xlab = "ISOMAP Dimension", ylab = "Correlation",
     main = "  ISOMAP Dimension  v.s. 
     Treatment Effect (unacceptable)",
     col = col.red.blue[trans(reord.cor)], pch = 16)
axis(1, at = 1:7, ord.dim.eff)
text(1:7, reord.cor, round(reord.cor, 2), cex = 1.2)
# the post[2] bigger, the prob. of unaccept. larger.
# so corr.(+), dim.val. big => unaccept.
#         (-)               =>   accept.
# legend
rasterImage(t(col.red.blue), 5, -0.9, 7, -0.8)
text(c(5, 6, 6, 7), c(-0.95, -0.7, -0.95, -0.95), 
     c(-1, "Correlation", 0, 1))

# b). dim.val.(isomap$points) => algorithm
#     Sig. algorithm
cor.dim.var <- apply(data.sig.var, 2, function(x)
  cor(x, lsvt.isomap$points))
# heatmap
heatmap(trans(t(cor.dim.var)), 
        col = col.red.blue,
        Rowv = NULL, 
        margins = c(4, 11),
        xlab = "ISOMAP Dimension", ylab = "Algorithm",
        main = "Algorithm v.s. ISOMAP Dimension")
# legend
rasterImage(t(col.red.blue), 5.4, 1, 6.7, 1.1)
text(c(5.28, 4.8, 6.05, 6.8), rep(1.05, 4), 
     c(-1, "Corr.", 0, 1))

# Medium or High Correlation is considered
positive.var <- apply(cor.dim.var, 1, function(x)
  names(which(x >= 0.3)))
negative.var <- apply(cor.dim.var, 1, function(x)
  names(which(x <= -0.3)))
# view~.
######################################################
# Dim1: 信號能量的穩定與否, 信號噪音比, 
#       小波波段1 & 3, 信號噪音比, 
#       聲帶振動頻率、振幅離散(偏離)趨勢
# + MFCC_1st,2nd,3rd, IMF, DFA
# - entropy_log_1~10, entropy_log3_1~5,8,9, HNR,
#   Jitter->F0/pitch, Shimmer->Ampl

# Dim2: 信號能量的穩定與否, 振幅偏離趨勢
#       信號能量的穩定與否, 信噪比, 小波波段3前段, 振幅
# + MFCC_0th, Log energy, Shimmer->Ampl_TKEO_prc25
# - MFCC_2nd,3rd, , IMF, DFA, entropy_log3_2~4
#   Shimmer->Ampl

# Dim3: 信號能量的穩定與否, 信噪比, 振幅
#       信噪比
# + MFCC_0th, IMF, Shimmer->Ampl_AM
# - VFER->SNR_SEO,SNR_TKEO1,NSR_TKEO1 

# Dim4: 信噪比
#       振幅偏離
# + VFER->NSR_TKEO1
# - Shimmer->Ampl_AM, abs0th_perturb

# Dim5: 信噪比
# + x
# - VFER->SNR_SEO,SNR_TKEO1,NSR_TKEO1

# Dim6: 
# + x
# - x

# Dim7: 
# + x
# - x
#######################################################
# 7. scores
iso.dim.def <- c("能量&音調持續穩定性",
                 "噪音比例", "音量波動性")
# a). 2D, risk (probability) of unaccept.
par(mfrow = c(2, 2))
par(mai = c(0.7, 0.7, 0.6, 0.2))
apply(matrix(c(1, 1, 2, 2, 3, 3), ncol = 2), 1, 
      function(x) {
        plot(lsvt.isomap$points[ , x[1]],
             lsvt.isomap$points[ , x[2]],
             col = col.red.blue[
               trans(lda.isomap$posterior[ , 2]*2-1)],
             pch = 16, cex = 1.2, 
             xlab = paste("Dim", x[1], iso.dim.def[x[1]]),
             ylab = paste("Dim", x[2], iso.dim.def[x[2]]), 
             cex.lab = 1.2)
        abline(h = 0, v = 0, lty = 2)
      }
)
plot(0, xlab = "", ylab = "", type = "n", xaxt = "n", 
     yaxt = "n", main = "\nProbability of Unacceptance",
     xlim = c(-1, 1), ylim = c(-1, 1))
rasterImage(t(col.red.blue), -0.5, -0.1, 0.5, 0.1)
text(c(-0.5, 0, 0.5, -0.5, 0.5),
     c(0.2, 0.2, 0.2, -0.2, -0.2), cex = 1.2,
     c(0, 0.5, 1, "Accept.", "Unaccept."))
title("\nISOMAP dimension projection", 
      outer = T, cex.main = 1.8)

# 3D
library(rgl)
apply(as.matrix((0:18)*20), 1, function(x) {
  #open3d()
  rgl.viewpoint(x, -10)
  plot3d(lsvt.isomap$points[ , 1:3], 
         main = "Isomap", 
         col = col.red.blue[
           trans(lda.isomap$posterior[ , 2]*2-1)], 
         size = 1.5, type = "s", xlab = "Dim.1",
         ylab = "Dim.2", zlab = "Dim.3")
  rgl.snapshot(paste("isomap_", x, ".png"), 
               fmt = "png", top = F)
})

# play3d(spin3d(axis = c(0, 0, 1), rpm = 10), 
#        duration = 10)


# b). 2D, Actual/pred. class
pred.class <- c(1 , 2)[lda.isomap$class]
actual.pred.class <- ifelse((class == 1 & 
                               pred.class == 2), 
                            3, ifelse((class == 2 & 
                                         pred.class == 1),
                                      4, class))
#   actual pred
# 1 accept accept
# 2 u      u
# 3 a      u
# 4 u      a

par(mfrow = c(2, 2))
par(mai = c(0.7, 0.7, 0.6, 0.2))
apply(matrix(c(1, 1, 2, 2, 3, 3), ncol = 2), 1, 
      function(x) {
        plot(lsvt.isomap$points[ , x[1]],
             lsvt.isomap$points[ , x[2]],
             col = c(4, 2, 1, 1)[actual.pred.class],
             pch = c(16, 17, 17, 16)[actual.pred.class],
             cex = 1.2, cex.lab = 1.2, 
             xlab = paste("Dim", x[1], iso.dim.def[x[1]]),
             ylab = paste("Dim", x[2], iso.dim.def[x[2]]))
        abline(h = 0, v = 0, lty = 2)
      }
)
plot(0, xlab = "", ylab = "", type = "n", xaxt = "n", 
     yaxt = "n", main = "\nActual/Prediction")
legend("center", c("Truly predicted Accept.",
                   "Truly predicted Unccept.",
                   "Falsely predicted Acccept.",
                   "Falsely predicted Unaccept."),
       col = c(4, 2, 1, 1), pch = c(16, 17, 17, 16),
       bty = "n", pt.cex = 1.2, cex = 1.2)
title("\nISOMAP dimension projection", 
      outer = T, cex.main = 1.8)
# 3D
apply(as.matrix((0:18)*20), 1, function(x) {
  open3d()
  rgl.viewpoint(x, -10)
  plot3d(lsvt.isomap$points[ , 1:3], 
         main = "Isomap", 
         col = c(4, 2, 1, 1)[actual.pred.class],
         size = 1.5, type = "s", xlab = "Dim.1",
         ylab = "Dim.2", zlab = "Dim.3")
  rgl.snapshot(paste("isomap_", x, ".png"), 
               fmt = "png", top = F)
})
rgl.snapshot("isomap_000.png", 
             fmt = "png", top = F)
# play3d(spin3d(axis = c(0, 0, 1), rpm = 10), 
#        duration = 10)


# 8. Unaccept. phonation sympton type
#　　unaccept. isomap$points clustered
# a). Unaccept. phonation isomap$points
unaccept.iso.pt <- lsvt.isomap$points[
  which(class == 2), ]

# b). Clustering Measurment
library(clValid)
intern <- clValid(unaccept.iso.pt, 2:10, 
                  clMethods = c("hierarchical", 
                                "kmeans", "pam"),
                  validation = "internal")
# summary(intern) # view results
stab <- clValid(unaccept.iso.pt, 2:10, 
                clMethods = c("hierarchical", 
                              "kmeans","pam"),
                validation = "stability") 
# summary(stab)

# c). Hierarchical, K = 2 
hclust.iso.pts <- hclust(dist(unaccept.iso.pt),
                         members = NULL)
# Clustered isomap pts plot, k = 2
cut.h.cluster <- cutree(hclust.iso.pts, k = 2)

par(mfrow = c(2, 2))
par(mai = c(0.7, 0.7, 0.6, 0.2))
apply(matrix(c(1, 1, 2, 2, 3, 3), ncol = 2), 1, 
      function(x) {
        plot(unaccept.iso.pt[ , x[1]],
             unaccept.iso.pt[ , x[2]],
             col = c("purple", "orange")[cut.h.cluster],
             pch = 16, cex = 1.2, 
             xlab = paste("Dim", x[1], iso.dim.def[x[1]]),
             ylab = paste("Dim", x[2], iso.dim.def[x[2]]), 
             cex.lab = 1.2)
        abline(h = 0, v = 0, lty = 2)
      }
)
plot(0, xlab = "", ylab = "", type = "n", 
     xaxt = "n", yaxt = "n", main = "",
     xlim = c(-1, 1), ylim = c(-1, 1))
legend("center", c("Group1 音調控制", 
                   "Group2 音量控制"),
       title = "未改善的症狀", pch = 16, 
       col = c("purple", "orange"), bty = "n",
       cex = 1.2, pt.cex = 1.2)
title("\nISOMAP dimension projection (Unacceptance)", 
      outer = T, cex.main = 1.8)


