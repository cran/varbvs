## ---- echo = FALSE, message = FALSE--------------------------------------
knitr::opts_chunk$set(eval = FALSE,collapse = TRUE,comment = "#")

## ---- eval = TRUE, message = FALSE---------------------------------------
library(lattice)
library(varbvs)

## ---- eval = TRUE--------------------------------------------------------
set.seed(1)

## ------------------------------------------------------------------------
#  load("cd.RData")

## ------------------------------------------------------------------------
#  r <- system.time(fit <- varbvs(X,NULL,y,"binomial",logodds = seq(-6,-3,0.25)))
#  cat(sprintf("Model fitting took %0.2f minutes.\n",r["elapsed"]/60))

## ------------------------------------------------------------------------
#  w   <- c(normalizelogweights(fit$logw))
#  pip <- c(varbvsindep(fit,X,NULL,y)$alpha %*% w)

## ------------------------------------------------------------------------
#  save(list = c("fit","map","pip","r"),
#       file = "varbvs.demo.cd.RData")

## ------------------------------------------------------------------------
#  print(summary(fit,nv = 9))

## ---- fig.width = 9,fig.height = 4,fig.align = "center"------------------
#  i <- which(fit$alpha %*% w > 0.5)
#  var.labels <- paste0(round(map$pos[i]/1e6,digits = 2),"Mb")
#  print(plot(fit,groups = map$chr,vars = i,var.labels = var.labels,gap = 7500,
#             ylab = "posterior prob."),
#        split = c(1,1,1,2),more = TRUE)
#  print(plot(fit,groups = map$chr,score = log10(pip + 0.001),vars = i,
#             var.labels = var.labels,gap = 7500,ylab = "log10 posterior prob."),
#        split = c(1,2,1,2),more = FALSE)

