library(jointCalib)
library(sampling)
library(laeken)
library(survey)
library(data.table)
library(ggplot2)
library(scales)
library(nonprobsvy)
library(xtable)

seed_for_sim <- 2023-07-14
set.seed(seed_for_sim+1)
n_reps <- 1000

## function taken from the publication Zhang, S., Han, P., & Wu, C. (2022). Calibration Techniques Encompassing Survey Sampling, Missing Data Analysis and Causal Inference. International Statistical Review.
cal.weight.EL <- function(g) {
  Fn <- function(rho, ghat) {
    -sum(log(1 + ghat %*% rho))
  }
  Grd <- function(rho, ghat) {
    -colSums(ghat / c(1 + ghat %*% rho))
  }
  rho.hat <-
    constrOptim(
      theta = rep(0, NCOL(g)),
      f = Fn,
      grad = Grd,
      ui = g,
      ci = rep(1 / NROW(g) - 1, NROW(g)),
      ghat = g
    )$par
  weight <- c(1 / (1 + g %*% rho.hat) / NROW(g))
  weight <- weight / sum(weight)
  return(weight)
}

# generate data according to the paper ------------------------------------

N <- 20000
n_r <- 10000
z1 <- rbinom(N, 1, 0.5)
z2 <- runif(N,0,2)
z3 <- rexp(N,1)
z4 <- rchisq(N, 4)
x1 <- z1
x2 <- z2 + 0.3*x1
x3 <- z3 + 0.2*(x1 + x2)
x4 <- z4 + 0.1*(x1 + x2 + x3)
e <- rnorm(N)

etol <- 1e-8
sigma03 <- stats::uniroot(f = function(s) cor(2+x1+x2+x3+x4, 2+x1+x2+x3+x4+s*e) - 0.3, c(0,20), tol = etol)$root
sigma05 <- stats::uniroot(f = function(s) cor(2+x1+x2+x3+x4, 2+x1+x2+x3+x4+s*e) - 0.5, c(0,20), tol = etol)$root
sigma08 <- stats::uniroot(f = function(s) cor(2+x1+x2+x3+x4, 2+x1+x2+x3+x4+s*e) - 0.8, c(0,20), tol = etol)$root

y1 <- 2+x1+x2+x3+x4+sigma03*e
y2 <- 2+x1+x2+x3+x4+sigma05*e
y3 <- 2+x1+x2+x3+x4+sigma08*e

theta_r <- stats::uniroot(f = function(s) sum(plogis(s + 0.1*x1+0.2*x2+0.3*x3+0.2*x4))-n_r, c(-20,0), tol = etol)$root
p_r <- plogis(theta_r + 0.1*x1+0.2*x2+0.3*x3+0.2*x4)

pop_data <- data.frame(x1,x2,x3,x4,y1,y2,y3,p_r) |> setDT()


# known population quantities ---------------------------------------------

x_totals <- with(pop_data, c("(Intercept)"=N, "x1"=sum(x1), "x2"=sum(x2), "x3"=sum(x3), "x4"=sum(x4)))
#p_quantiles <- seq(0.25, 0.75, 0.25) ## for estimation
p_quantiles <- c(seq(0.25, 0.75, 0.25), seq(0.1, 0.9, 0.1))
#p_quantiles <- seq(0.1, 0.9, 0.1) 
#p_quantiles <- p_quantiles_est ## for calibration
x2_q <- with(pop_data, quantile(x2, p_quantiles))
x3_q <- with(pop_data, quantile(x3, p_quantiles))
x4_q <- with(pop_data, quantile(x4, p_quantiles))
x_quants <- list(x2=x2_q,x3=x3_q, x4=x4_q)

# known population y ------------------------------------------------------

pop_true_vals <- pop_data[, lapply(.SD, 
                                   function(x) data.frame(mean=mean(x), 
                                                          q25=quantile(x, 0.25), 
                                                          q50=quantile(x, 0.5), 
                                                          q75=quantile(x, 0.75))), 
                          .SDcols = c("y1", "y2", "y3")] |> 
  melt(value.name = "true") |> 
  {\(x) x[, c("y", "stat") := tstrsplit(variable, "\\.")][, .(y, stat, true)]}()


# simulation --------------------------------------------------------------
results <- list()

for (r in 1:n_reps) {
  set.seed(r)
  if (r %% 50 == 0) cat("Iteration: ", r, "\n")
  
  ## nonprob sample
  sample_resp <- pop_data[which(sampling::UPpoisson(p_r)==1),]
  sample_resp[, w_naive:= N/.N]
  
  ## IPW
  ipw <- nonprob(selection = ~ x1 + x2 + x3 + x4,
                 target = ~ y1,
                 data = sample_resp,
                 pop_totals = x_totals,
                 control_selection = controlSel(est_method_sel = "gee", h_x = "2"))
  
  ## calibration with quantiles only
  w_res_d <- joint_calib(formula_quantiles = ~ x2 + x3 + x4,
                         data = sample_resp,
                         dweights = sample_resp$w_naive,
                         N = N,
                         pop_quantiles = x_quants,
                         method = "raking",
                         backend = "sampling",
                         control = control_calib(interpolation = "linear") # linear
  )
  
  w_res <- w_res_d$g
  
  
  ## calibration with quantiles and totals
  w_res2_d <- joint_calib(formula_quantiles = ~ x2 + x3 + x4,
                          formula_totals = ~ x1 + x2 + x3 + x4,
                          data = sample_resp,
                          dweights = sample_resp$w_naive,
                          N = N,
                          pop_quantiles = x_quants,
                          pop_totals = x_totals[-1],
                          method = "raking",
                          backend = "sampling",
                          control = control_calib(interpolation = "linear") # linear
  )
  
  w_res2 <- w_res2_d$g
  
  ## calibration with totals only
  w_res_st <- calib(Xs = with(sample_resp, cbind(1, x1, x2, x3, x4)), 
                    d = sample_resp$w_naive, 
                    total = x_totals, 
                    method = "raking")
  
  ## el with means only
  A_diff0 <- w_res2_d$Xs[,c("x1", "x2", "x3", "x4")] - 
    do.call("rbind", replicate(nrow(w_res2_d$Xs), w_res2_d$totals[c(c("x1", "x2", "x3", "x4"))]/N, simplify = F))
  w_res_el0 <- cal.weight.EL(g = A_diff0)
  
  ## quantiles only
  A_diff1 <- w_res_d$Xs[,-1] - do.call("rbind", replicate(nrow(w_res_d$Xs), w_res_d$totals[-1]/N, simplify = F))
  w_res_el1 <- cal.weight.EL(g = A_diff1)
  
  ## quantiles and means
  A_diff2 <- w_res2_d$Xs[,-1] - do.call("rbind", replicate(nrow(w_res2_d$Xs), w_res2_d$totals[-1]/N, simplify = F))
  w_res_el2 <- cal.weight.EL(g = A_diff2)
  
  ## means only
  
  ## weights
  sample_resp[, ":="(wcal0=w_naive*w_res_st, qcal1=w_naive*w_res, qcal2=w_naive*w_res2, wipw=ipw$weights, 
                     wel0=w_res_el0*N, wel1=w_res_el1*N, wel2=w_res_el2*N)]
  ## means
  results_mean_n <- sample_resp[, lapply(.SD, mean), .SDcols = patterns("y")][, w:="naive"][, stat:="mean"]
  results_mean_wcal0 <- sample_resp[, lapply(.SD, weighted.mean, w = wcal0), .SDcols = patterns("y")][, w:="wcal0"][, stat:="mean"]
  results_mean_qcal1 <- sample_resp[, lapply(.SD, weighted.mean, w = qcal1), .SDcols = patterns("y")][, w:="qcal1"][, stat:="mean"]
  results_mean_qcal2 <- sample_resp[, lapply(.SD, weighted.mean, w = qcal2), .SDcols = patterns("y")][, w:="qcal2"][, stat:="mean"]
  results_mean_ipw <- sample_resp[, lapply(.SD, weighted.mean, w = wipw), .SDcols = patterns("y")][, w:="wipw"][, stat:="mean"]
  results_mean_el0 <- sample_resp[, lapply(.SD, weighted.mean, w = wel0), .SDcols = patterns("y")][, w:="wel0"][, stat:="mean"]
  results_mean_el1 <- sample_resp[, lapply(.SD, weighted.mean, w = wel1), .SDcols = patterns("y")][, w:="wel1"][, stat:="mean"]
  results_mean_el2 <- sample_resp[, lapply(.SD, weighted.mean, w = wel2), .SDcols = patterns("y")][, w:="wel2"][, stat:="mean"]
  ## q25
  results_q25_n <- sample_resp[, lapply(.SD, quantile, probs=0.25), .SDcols = patterns("y")][, w:="naive"][, stat:="q25"]
  results_q25_wcal0 <- sample_resp[, lapply(.SD, weightedQuantile, probs=0.25, w = wcal0), .SDcols = patterns("y")][, w:="wcal0"][, stat:="q25"]
  results_q25_qcal1 <- sample_resp[, lapply(.SD, weightedQuantile, probs=0.25, w = qcal1), .SDcols = patterns("y")][, w:="qcal1"][, stat:="q25"]
  results_q25_qcal2 <- sample_resp[, lapply(.SD, weightedQuantile, probs=0.25, w = qcal2), .SDcols = patterns("y")][, w:="qcal2"][, stat:="q25"]
  results_q25_ipw <-  sample_resp[, lapply(.SD, weightedQuantile, probs=0.25, w = wipw), .SDcols = patterns("y")][, w:="wipw"][, stat:="q25"]
  results_q25_el0 <-  sample_resp[, lapply(.SD, weightedQuantile, probs=0.25, w = wel0), .SDcols = patterns("y")][, w:="wel0"][, stat:="q25"]
  results_q25_el1 <-  sample_resp[, lapply(.SD, weightedQuantile, probs=0.25, w = wel1), .SDcols = patterns("y")][, w:="wel1"][, stat:="q25"]
  results_q25_el2 <-  sample_resp[, lapply(.SD, weightedQuantile, probs=0.25, w = wel2), .SDcols = patterns("y")][, w:="wel2"][, stat:="q25"]

  ## q50
  results_q50_n <- sample_resp[, lapply(.SD, quantile, probs=0.50), .SDcols = patterns("y")][, w:="naive"][, stat:="q50"]
  results_q50_wcal0 <- sample_resp[, lapply(.SD, weightedQuantile, probs=0.5, w = wcal0), .SDcols = patterns("y")][, w:="wcal0"][, stat:="q50"]
  results_q50_qcal1 <- sample_resp[, lapply(.SD, weightedQuantile, probs=0.5, w = qcal1), .SDcols = patterns("y")][, w:="qcal1"][, stat:="q50"]
  results_q50_qcal2 <- sample_resp[, lapply(.SD, weightedQuantile, probs=0.5, w = qcal2), .SDcols = patterns("y")][, w:="qcal2"][, stat:="q50"]
  results_q50_ipw <-  sample_resp[, lapply(.SD, weightedQuantile, probs=0.5, w = wipw), .SDcols = patterns("y")][, w:="wipw"][, stat:="q50"]
  results_q50_el0 <-  sample_resp[, lapply(.SD, weightedQuantile, probs=0.5, w = wel0), .SDcols = patterns("y")][, w:="wel0"][, stat:="q50"]
  results_q50_el1 <-  sample_resp[, lapply(.SD, weightedQuantile, probs=0.5, w = wel1), .SDcols = patterns("y")][, w:="wel1"][, stat:="q50"]
  results_q50_el2 <-  sample_resp[, lapply(.SD, weightedQuantile, probs=0.5, w = wel2), .SDcols = patterns("y")][, w:="wel2"][, stat:="q50"]
  ## q75
  results_q75_n <- sample_resp[, lapply(.SD, quantile, probs=0.75), .SDcols = patterns("y")][, w:="naive"][, stat:="q75"]
  results_q75_wcal0 <- sample_resp[, lapply(.SD, weightedQuantile, probs=0.75, w = wcal0), .SDcols = patterns("y")][, w:="wcal0"][, stat:="q75"]
  results_q75_qcal1 <- sample_resp[, lapply(.SD, weightedQuantile, probs=0.75, w = qcal1), .SDcols = patterns("y")][, w:="qcal1"][, stat:="q75"]
  results_q75_qcal2 <- sample_resp[, lapply(.SD, weightedQuantile, probs=0.75, w = qcal2), .SDcols = patterns("y")][, w:="qcal2"][, stat:="q75"]
  results_q75_ipw <-  sample_resp[, lapply(.SD, weightedQuantile, probs=0.75, w = wipw), .SDcols = patterns("y")][, w:="wipw"][, stat:="q75"]
  results_q75_el0 <-  sample_resp[, lapply(.SD, weightedQuantile, probs=0.75, w = wel0), .SDcols = patterns("y")][, w:="wel0"][, stat:="q75"]
  results_q75_el1 <-  sample_resp[, lapply(.SD, weightedQuantile, probs=0.75, w = wel1), .SDcols = patterns("y")][, w:="wel1"][, stat:="q75"]
  results_q75_el2 <-  sample_resp[, lapply(.SD, weightedQuantile, probs=0.75, w = wel2), .SDcols = patterns("y")][, w:="wel2"][, stat:="q75"]
  
  results[[r]] <- rbind(results_mean_n, results_mean_wcal0, results_mean_qcal1, results_mean_qcal2, results_mean_ipw, 
                        results_mean_el0, results_mean_el1, results_mean_el2,
                        results_q25_n, results_q25_wcal0, results_q25_qcal1, results_q25_qcal2, results_q25_ipw,
                        results_q25_el0,results_q25_el1,results_q25_el2,
                        results_q50_n, results_q50_wcal0, results_q50_qcal1, results_q50_qcal2, results_q50_ipw,
                        results_q50_el0, results_q50_el1,results_q50_el2,
                        results_q75_n, results_q75_wcal0, results_q75_qcal1, results_q75_qcal2, results_q75_ipw,
                        results_q75_el0,results_q75_el1,results_q75_el2)
}


results_df <- rbindlist(results, idcol = "r") |> 
  melt(id.vars = c("r", "w", "stat"), variable.name = "y") |> 
  transform(y=as.character(y))

results_df[pop_true_vals, on = c("y", "stat"), true:=i.true]
results_df[, w:=factor(w, 
                       c("naive", "wipw", "wcal0", "qcal1", "qcal2", "wel0", "wel1","wel2"), 
                       c("Naive", "IPW", "CAL", "QCAL1", "QCAL2", "EL", "QEL1","QEL2"))]

saveRDS(results_df, file = "results/sim-chen2020-jasa-modif.rds")
fwrite(x = results_df, file = "results/sim-chen2020-jasa-modif.csv.gz")

