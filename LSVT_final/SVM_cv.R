## Feature & Dimension Reduction method Selection
# (criterion: Max. classification(SVM) accuracy)

# 1. Data updated by 
#    Oredered significant variables (algorithm)

# d.algor.ord <- as.matrix(data)[ , order.bw]

# 2. SVM, CV(Leave one out), function
library(e1071)
cv.svm.accur.fun <- function(test.id, data, dim) {
  x.train <- data[-test.id, 1:dim]
  x.test <- data[test.id, 1:dim]
  y.train <- class.f[-test.id]
  y.test <- class.f[test.id]
  train.set <- data.frame(x.train, y.train)
  
  svm.model <- svm(y.train ~ ., data = train.set)
  
  pred <- predict(svm.model, x.test)
  confus.mat <- table(pred, y.test)
  accuracy <- sum(diag(prop.table(confus.mat)))
  accuracy
}

# 2. Dimension Reduction method Selection (matplot.)
library(dplyr)
cv.algor.DR.method.select.svm <- function(algor.num) {
  LSVT <- d.algor.ord[ , 1:algor.num]
  # a). PC scores
  # library(FactoMineR)
  # library(ggplot2)
  # library(factoextra)
  lsvt.pca <- PCA(LSVT, scale.unit = T, graph = F, 
                  ncp = 10)
  eig <- get_eigenvalue(lsvt.pca)
  var <- get_pca_var(lsvt.pca)
  ind <- get_pca_ind(lsvt.pca)
  lsvt.pc.score <- ind$coord[ , 1:10]
  svm.pc.accur <- apply(as.matrix(2:10), 1, function(x)
    apply(as.matrix(1:dim(LSVT)[1]), 1, cv.svm.accur.fun,
          data.frame(lsvt.pc.score), x) %>% mean)
  
  # b). MDS projected value
  mds.subj <- cmdscale(dist(scale(LSVT)), k = 10)
  svm.mds.accur <- apply(as.matrix(2:10), 1, function(x)
    apply(as.matrix(1:dim(LSVT)[1]), 1, cv.svm.accur.fun,
          data.frame(mds.subj), x) %>% mean)
  
  # c). ISOMAP projected value
  # library(vegan)
  iso.score.k <- apply(as.matrix(2:20), 1, function(neighbor) {
    lsvt.isomap <- isomap(dist(scale(LSVT)), 
                          ndim = 10, k = neighbor)$points
    dim.accur <- apply(as.matrix(2:10), 1, function(x)
      apply(as.matrix(1:dim(LSVT)[1]), 1, cv.svm.accur.fun,
            data.frame(lsvt.isomap), x) %>% mean)
    dim.accur
  })
  svm.iso.accur <- apply(iso.score.k, 1, max)
  
  svm.accur <- rbind(rep(NA, 3),
                     data.frame(PCA = svm.pc.accur,
                                MDS = svm.mds.accur,
                                ISOMAP = svm.iso.accur))
  matplot(1 - svm.accur, pch = c(1, 3, 2), type = "b",
          col = 2:4, lwd = 2, cex = c(2, 1.5, 1.2),
          main = paste("\n", algor.num, " algorithms"),
          xlab = "Dimensions", ylab = "Error rate")
}
# Before feature selection
cv.algor.DR.method.select.svm(dim(d.algor.ord)[2])
legend("topright", c("PCA", "MDS", "ISOMAP"),
       pch = c(1, 3, 2), lwd = 2, col = 2:4,
       lty = 2:4, pt.cex = c(2, 1.5, 1.2))
# After feature selection (ISOMAP is always the best)
par(mfrow = c(3, 3))
par(mai = c(0.6, 0.6, 0.5, 0.2))
cv.svm.algor.DR.method.select <- 
  apply(as.matrix(c(20, 40, 60, 80, 100, 150, 200, 250, 309)),
        1, cv.algor.DR.method.select.svm)
legend("top", c("PCA", "MDS", "ISOMAP"),
       pch = c(1, 3, 2), lwd = 2, col = 2:4, horiz = T,
       lty = 1:3, pt.cex = c(1.5, 1.2, 1.2), cex = 0.8,
       bg = cm.colors(1, alpha = 0), bty = "n")
title("\nError Rate by SVM(CV) upon Variables Selection",
      outer = T)

# 3. Feature selection
# Compare isomap error rate in diff. number of var. selected
cv.svm.algor.isomap.select.fun <- function(algor.num) {
  LSVT <- d.algor.ord[ , 1:algor.num]
  # c). ISOMAP projected value
  iso.score.k <- apply(as.matrix(2:20), 1, function(neighbor) {
    lsvt.isomap <- isomap(dist(scale(LSVT)), 
                          ndim = 10, k = neighbor)$points
    dim.accur <- apply(as.matrix(2:10), 1, function(x)
      apply(as.matrix(1:dim(LSVT)[1]), 1, cv.svm.accur.fun,
            data.frame(lsvt.isomap), x) %>% mean)
    dim.accur
  })
  svm.iso.accur <- apply(iso.score.k, 1, max)
  svm.iso.accur
}
cv.svm.algor.isomap.select <- 
  rbind(rep(NA, 9),
        apply(as.matrix(c(20, 40, 60, 80, 100, 
                          150, 200, 250, 309)),
              1, cv.svm.algor.isomap.select.fun))

library(fields)
matplot(1-cv.svm.algor.isomap.select, pch = 16,
        type = "l", lwd = 4, lty = 1,
        xlim = c(0.5, 10.5), ylim = c(0.06, 0.28),
        col = tim.colors(9), xlab = "Dimensions",
        ylab = "Error rate", 
        main = "ISOMAP Classification Error Rate 
        (Leave-one-out SVM)",
        xaxt = "n")
axis(1, at = 2:10, 2:10)
legend("topright", legend = rev(c(20, 40, 60, 80, 100, 150, 
                                  200, 250, 309)),
       lwd = 3, col = rev(tim.colors(9)), ncol = 3,
       title = "Feature selection", cex = 0.8)
points(rep(6, 2), 
       rep(1-cv.svm.algor.isomap.select[6, 2], 2),
       pch = c(16, 1), col = 2, cex = c(1, 3)) # (Continue...)
text(4.7, 1-cv.svm.algor.isomap.select[6, 2]-0.01,
     paste("Accuracy\n= ", round(100*cv.svm.algor.isomap.select[6, 2]), "%"),
     col = 2, cex = 1.2)
# Thus, the first 40 sig. var.(features) are selected
# to classify trt. consequences (SVM) 
# after Dimension Reduction (ISOMAP dim.6).

# 4. Data that contains only the 40 most sig. var.
# data.sig.var <- d.algor.ord[ , 1:40]

# 5. SVM, CV(Leave one out), isomap dim=6
# a). Find the most proper k(# of neighbors)
library(dplyr)
cv.svm.accur.given.k <- apply(as.matrix(2:20), 1, function(x) {
  iso.pts <- isomap(dist(scale(data.sig.var)), 
                    ndim = 6, k = x)$points
  accuracy <- apply(as.matrix(1:dim(data.sig.var)[1]), 1,
                    cv.svm.accur.fun,
                    data.frame(iso.pts), 6) %>% mean
  accuracy 
})
cv.svm.best.k <- (2:20)[which.max(cv.svm.accur.given.k)]

# b). Isomap
cv.svm.lsvt.isomap <- isomap(dist(scale(data.sig.var)), 
                             ndim = 6, k = cv.svm.best.k)
# c). classification by SVM
library(e1071)
cv.svm.prob.fun <- function(test.id) {
  data <- data.frame(cv.svm.lsvt.isomap$points)
  x.train <- data[-test.id, ]
  x.test <- data[test.id, ]
  y.train <- class.f[-test.id]
  y.test <- class.f[test.id]
  train.set <- data.frame(x.train, y.train)
  
  svm.model <- svm(y.train ~ ., data = train.set,
                   probability = T)
  
  pred <- predict(svm.model, x.test, probability = T)
  c(attr(pred, "probabilities"), pred)
}
cv.svm.prob.fun.out <- t(apply(as.matrix(1:dim(data.sig.var)[1]), 
                               1, cv.svm.prob.fun))
cv.svm.prob.class <- apply(cv.svm.prob.fun.out, 1, function(x) {
  if ((x[3] == 1) && (x[1] < x[2]) || 
      (x[3] == 2) && (x[1] > x[2])) {
    replace <- x[1]
    x[1] <- x[2]
    x[2] <- replace
  }
  x
}) %>% t
colnames(cv.svm.prob.class) <- c("Acceptable",
                                 "Unacceptable",
                                 "class")
# 6. Prob.of unaccept. [corr.to] iso.dim.
cor.dim.unaccept.prob <-
  apply(cv.svm.lsvt.isomap$points, 2, function(dim)
    cor(dim, cv.svm.prob.class[ , 2]))

plot(1:6, cor.dim.unaccept.prob, 
     cex = 20 * (abs(cor.dim.unaccept.prob)),
     xlim = c(0.5, 6.5), ylim = c(-1, 1),
     xlab = "ISOMAP Dimension", ylab = "Correlation",
     main = "  ISOMAP Dimension 
     v.s.  Unacceptance in Probability",
     col = col.red.blue[trans(cor.dim.unaccept.prob)], pch = 16)
abline(h = 0, lty = 2)
text(1:6, cor.dim.unaccept.prob, 
     round(cor.dim.unaccept.prob, 2), cex = 1.2)
rasterImage(t(col.red.blue), 4, -0.9, 6, -0.8)
text(c(4, 5, 5, 6), c(-0.95, -0.7, -0.95, -0.95), 
     c(-1, "Correlation", 0, 1))

# 7. isomap$points [corr.to] algorithm
#     Sig. algorithm
cv.svm.iso.dim <- cv.svm.lsvt.isomap$points

library(fields)
col.red.blue <- two.colors(start = "blue", 
                           middle = "gray", 
                           end = "red")
trans <- function(x) {round((x+1)*127.5)+1}

cor.dim.var <- apply(data.sig.var, 2, function(x)
  cor(x, cv.svm.iso.dim))
# heatmap
heatmap(trans(t(cor.dim.var)), 
        col = col.red.blue,
        Rowv = NULL, 
        margins = c(4, 11),
        xlab = "ISOMAP Dimension", ylab = "Algorithm",
        main = "Algorithm v.s. ISOMAP Dimension")
# legend
rasterImage(t(col.red.blue), 5.1, 1, 6.4, 1.1)
text(c(4.98, 4.5, 5.7, 6.5), rep(1.05, 4), 
     c(-1, "Corr.", 0, 1))

# High Correlation is considered
cv.svm.positive.var <- apply(cor.dim.var, 1, function(x)
  names(which(x >= 0.5)))
cv.svm.negative.var <- apply(cor.dim.var, 1, function(x)
  names(which(x <= -0.5)))
# view~. => def. dim.
iso.dim.def <- c("能量、音調持續穩定性",
                 "噪音(聲帶)、低音量波動性", 
                 "噪音(聲帶、咬合)")

# 8. plot
# Density plot (iso.pts' distr. given pred.class)
par(mfrow = c(2, 3))
par(mai = c(0.3, 0.6, 0.6, 0.2))
apply(as.matrix(1:6), 1, function(x) {
  dim <- cv.svm.lsvt.isomap$points[ , x]
  # Actual 
  sm.density.compare(dim, class.f,
                     col = c("green", "yellow")[class],
                     lty = rep(1, 2), lwd = 4,
                     xlim = c(-30, 20), 
                     ylim = c(0, 0.2), xaxt = "n",
                     yaxt = "n", bty = "n",
                     xlab = "", ylab = "")
  par(new = T)
  # Pred.
  sm.density.compare(dim, cv.svm.prob.class[ , 3],
                     col = c(4, 2)[cv.svm.prob.class[ , 3]],
                     lty = rep(2, 2), lwd = 3,
                     xlim = c(-30, 20), 
                     ylim = c(0, 0.2),
                     xlab = "")
  titles <- ifelse(is.na(iso.dim.def[x]) == 1,
                   paste(paste("\n\nDim", x)), 
                   paste("\n\nDim", x, iso.dim.def[x]))
  title(titles)
})

legend("topleft", c("Accept./Pred.", "Accept./Actual",
                    "Unaccept./Pred.", "Unaccept./Actual"),
       col = c("blue", "green", "red", "yellow"),
       lty = c(2, 1, 2, 1), lwd = c(3, 4, 3, 4),
       bg = cm.colors(1, alpha = 0), cex = 1.1)
title("\nISOMAP dimensions density plot",
      outer = T, cex.main = 2)

# 2D prob.
par(mfrow = c(2, 2))
par(mai = c(0.7, 0.7, 0.6, 0.2))
apply(matrix(c(1, 1, 2, 2, 3, 3), ncol = 2), 1, 
      function(x) {
        plot(cv.svm.iso.dim[ , x[1]],
             cv.svm.iso.dim[ , x[2]],
             col = col.red.blue[
               trans(cv.svm.prob.class[ , 2]*2-1)],
             pch = 16, cex = 1.2, 
             xlab = paste("Dim", x[1], iso.dim.def[x[1]]),
             ylab = paste("Dim", x[2], iso.dim.def[x[2]]), 
             cex.lab = 1.1)
        abline(h = 0, v = 0, lty = 2)
      }
)
plot(0, xlab = "", ylab = "", type = "n", xaxt = "n", 
     yaxt = "n", main = "\nProbability of Unacceptance",
     xlim = c(-1, 1), ylim = c(-1, 1))
rasterImage(t(col.red.blue), -0.5, -0.1, 0.5, 0.1)
text(c(-0.5, 0, 0.5, -0.5, 0.5),
     c(0.22, 0.22, 0.22, -0.22, -0.22), cex = 1.2,
     c(0, 0.5, 1, "Accept.", "Unaccept."))
title("\nISOMAP dimension projection", 
      outer = T, cex.main = 1.7)

# 3D
library(rgl)
plot3d(cv.svm.iso.dim[ , 1:3], 
       type = "s", col = col.red.blue[
         trans(cv.svm.prob.class[ , 2]*2-1)],
       size = 1.3, xlab = "Dim.1",
       ylab = "Dim.2", zlab = "Dim.3")

# 2D, pred./actual class
prop.table(table(class, pred.class), 2)
pred.class <- cv.svm.prob.class[ , 3]
cv.svm.actual.pred.class <- 
  ifelse((class == 1 & pred.class == 2), 3, 
         ifelse((class == 2 & pred.class == 1), 4, class))

par(mfrow = c(2, 2))
par(mai = c(0.7, 0.7, 0.6, 0.2))
apply(matrix(c(1, 1, 2, 2, 3, 3), ncol = 2), 1, 
      function(x) {
        plot(cv.svm.iso.dim[ , x[1]],
             cv.svm.iso.dim[ , x[2]],
             col = c(4, 2, 1, 1)[cv.svm.actual.pred.class],
             pch = c(16, 17, 17, 16)[cv.svm.actual.pred.class],
             cex = 1.2, cex.lab = 1.1,
             xlab = paste("Dim", x[1], iso.dim.def[x[1]]),
             ylab = paste("Dim", x[2], iso.dim.def[x[2]]))
        abline(h = 0, v = 0, lty = 2)
      }
)
plot(c(-0.5, 0.5, -0.5, 0.5),
     c(0.5, 0.5, -0.5, -0.5),
     pch = c(16, 17, 16, 17),
     col = c(4, 1, 1, 2),
     xlab = "", ylab = "", cex = 2,
     xaxt = "n", yaxt = "n", main = "", 
     xlim = c(-1, 1), ylim = c(-1, 1))
abline(v = 0, h = 0)
mtext(c("Actual", "Prediction"), 
      side = 2:3, padj = c(-1.5, -1.5),
      adj = c(0.5, 0.5), cex = 1.2)
mtext(c("Accept.", "Unaccept.",
        "Acceptable", "Unacceptable"), 
      side = c(2 , 2, 3, 3),
      padj = rep(-0.2, 4), 
      adj = c(1, 0, 0.1, 0.9))
title("\nISOMAP dimension projection", 
      outer = T, cex.main = 1.8)

# 3D, accept.
plot3d(cv.svm.iso.dim[which(pred.class == 1), 1:3], 
       type = "s", xlab = "Dim.1",
       ylab = "Dim.2", zlab = "Dim.3",
       col = c(4, 2, 1, 1)[cv.svm.actual.pred.class][
         which(pred.class == 1)], size = 2)
# 3D, unaccept.
plot3d(cv.svm.iso.dim[which(pred.class == 2), 1:3], 
       type = "s", xlab = "Dim.1",
       ylab = "Dim.2", zlab = "Dim.3",
       col = c(4, 2, 1, 1)[cv.svm.actual.pred.class][
         which(pred.class == 2)], size = 2)


# 9. isomap$points clustered
# Clustered isomap pts plot, k = 4
hclust.iso.pts <- hclust(dist(cv.svm.iso.dim),
                         members = NULL)
cut.h.cluster <- cutree(hclust.iso.pts, k = 4)

par(mfrow = c(2, 2))
par(mai = c(0.7, 0.7, 0.6, 0.2))
apply(matrix(c(1, 1, 2, 2, 3, 3), ncol = 2), 1, 
      function(x) {
        plot(cv.svm.iso.dim[ , x[1]],
             cv.svm.iso.dim[ , x[2]],
             col = c("purple", "orange", 
                     "green", "red")[cut.h.cluster],
             pch = 16, cex = 1.2, 
             xlab = paste("Dim", x[1], iso.dim.def[x[1]]),
             ylab = paste("Dim", x[2], iso.dim.def[x[2]]), 
             cex.lab = 1)
        abline(h = 0, v = 0, lty = 2)
      }
)
plot(0, xlab = "", ylab = "", type = "n", 
     xaxt = "n", yaxt = "n", main = "",
     xlim = c(-1, 1), ylim = c(-1, 1))
legend("center", c("Group1 有改善", 
                   "Group2 症狀輕微",
                   "Group3 症狀嚴重"),
       pch = 16, bty = "n", cex = 1.2, pt.cex = 1.2, 
       col = c("purple", "orange", "green"))
title("\nISOMAP dimension projection", 
      outer = T, cex.main = 1.7)