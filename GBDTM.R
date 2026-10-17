
library(readxl)
library(trajeR)

#单轨迹
data <- as.data.frame(read_excel(file.choose(), sheet = 1))
indicator <- "DQ"                  
order <- c(2, 2, 2)               
limits <- c(0, 168)                

Y <- as.matrix(data[paste0(indicator, 1:4)])
keep <- rowSums(!is.na(Y)) > 0
Y <- Y[keep, , drop = FALSE]
ID <- data$ID[keep]
stopifnot(is.numeric(Y), !anyNA(ID), !anyDuplicated(ID),
          length(order) >= 2, all(order %in% 0:3))
stopifnot(all(Y >= limits[1] & Y <= limits[2], na.rm = TRUE))

months <- c(12, 18, 24, 36)
A <- matrix(rep((months - 24) / 12, each = nrow(Y)), nrow(Y))

fit <- trajeR(Y = Y, A = A, degre = order, Model = "CNORM",
              Method = "L", ssigma = TRUE,
              ymin = limits[1], ymax = limits[2],
              hessian = FALSE, itermax = 3000)
if (!fit$convergence %in% c(1, 2)) warning("模型未可靠收敛，请检查拟合结果。")

posterior <- GroupProb(fit, Y = Y, A = A)
stopifnot(all(is.finite(posterior)), max(abs(rowSums(posterior) - 1)) < 1e-6)
individual_results <- data.frame(ID, group = max.col(posterior, ties.method = "first"), posterior)
pi <- exp(fit$theta - max(fit$theta)); pi <- pi / sum(pi)
group_results <- data.frame(
  group = seq_along(order), proportion = pi,
  n = tabulate(individual_results$group, nbins = length(order)),
  AvePP = sapply(seq_along(order), function(g)
    mean(posterior[individual_results$group == g, g]))
)


n_parameters <- sum(order + 1) + 1 + length(order) - 1
model_results <- data.frame(logLik = fit$Likelihood,
  BIC = -2 * fit$Likelihood + n_parameters * log(nrow(Y)))
print(model_results)
print(group_results)
print(head(individual_results, 10))

colors <- hcl.colors(length(order), "Dark 3")
plotrajeR(fit, col = rep(colors, 2), xaxt = "n",
          xlab = "Age (months)", ylab = indicator)
axis(1, at = (months - 24) / 12, labels = months)
legend("topright", legend = paste0("G", seq_along(order)), col = colors, lwd = 2, bty = "n")

#双轨迹
library(readxl)
library(trajeR)
stopifnot(packageVersion("trajeR") >= "1.0")

data <- as.data.frame(read_excel(file.choose(), sheet = 1))
indicators <- c("DQ", "ZBMI")       
order1 <- c(2, 2, 2)              
order2 <- c(2, 2, 2)              
limits1 <- c(0, 168)               
limits2 <- c(-4, 4)

Y1 <- as.matrix(data[paste0(indicators[1], 1:4)])
Y2 <- as.matrix(data[paste0(indicators[2], 1:4)])
keep <- rowSums(!is.na(Y1)) + rowSums(!is.na(Y2)) > 0
Y1 <- Y1[keep, , drop = FALSE]
Y2 <- Y2[keep, , drop = FALSE]
ID <- data$ID[keep]
stopifnot(is.numeric(Y1), is.numeric(Y2), !anyNA(ID), !anyDuplicated(ID),
          length(order1) >= 2, length(order2) >= 2, all(c(order1, order2) %in% 0:3))
stopifnot(all(Y1 >= limits1[1] & Y1 <= limits1[2], na.rm = TRUE),
          all(Y2 >= limits2[1] & Y2 <= limits2[2], na.rm = TRUE))

months <- c(12, 18, 24, 36)
A <- matrix(rep((months - 24) / 12, each = nrow(Y1)), nrow(Y1))
make_formula <- function(y, degree) {
  rhs <- if (degree == 0) "1" else
    paste(c("A", if (degree >= 2) paste0("I(A^", 2:degree, ")")), collapse = " + ")
  as.formula(paste(y, "~", rhs))
}
formulas <- list(lapply(order1, function(d) make_formula("Y1", d)),
                 lapply(order2, function(d) make_formula("Y2", d)))

fit <- trajeR_mult(formulas = formulas,
                   ng = list(length(order1), length(order2)), models = c("CNORM", "CNORM"),
                   ssigma = TRUE, ymin = list(limits1[1], limits2[1]),
                   ymax = list(limits1[2], limits2[2]), hessian = FALSE,
                   control = list(maxit = 3000, reltol = 1e-10))
if (fit$convergence != 0) warning("模型未收敛，请检查拟合结果。")

probabilities <- get_group_probabilities(fit,
                                         lY = list(Y1, Y2), lA = list(A, A),
                                         lymin = list(limits1[1], limits2[1]), lymax = list(limits1[2], limits2[2]))
posterior <- probabilities$post_prob
stopifnot(all(is.finite(posterior)), max(abs(rowSums(posterior) - 1)) < 1e-6)
g1 <- length(order1); g2 <- length(order2)
grid <- expand.grid(first = seq_len(g1), second = seq_len(g2))
post1 <- sapply(seq_len(g1), function(g) rowSums(posterior[, grid$first == g, drop = FALSE]))
post2 <- sapply(seq_len(g2), function(g) rowSums(posterior[, grid$second == g, drop = FALSE]))
individual_results <- data.frame(ID,
                                 group1 = max.col(post1, ties.method = "first"), group2 = max.col(post2, ties.method = "first"),
                                 max_probability1 = apply(post1, 1, max), max_probability2 = apply(post2, 1, max))
names(individual_results)[2:3] <- paste0(indicators, "_group")

joint_probability <- matrix(probabilities$prior_prob[1, ], nrow = g1, ncol = g2,
                            dimnames = list(paste0(indicators[1], "_G", seq_len(g1)), paste0(indicators[2], "_G", seq_len(g2))))
P_second_given_first <- sweep(joint_probability, 1, rowSums(joint_probability), `/`)
P_first_given_second <- sweep(joint_probability, 2, colSums(joint_probability), `/`)


n_parameters <- sum(order1 + 1) + sum(order2 + 1) + 2 + g1 * g2 - 1
model_results <- data.frame(logLik = fit$likelihood,
                            BIC = -2 * fit$likelihood + n_parameters * log(nrow(Y1)))
print(model_results)
print(joint_probability)
print(P_second_given_first)         # 每行合计为 1：给定第一个指标组，第二个指标的组概率
print(P_first_given_second)         # 每列合计为 1：给定第二个指标组，第一个指标的组概率
print(head(individual_results, 10))

old_par <- par(mfrow = c(1, 2))
age_grid <- seq(12, 36, length.out = 100)
for (m in 1:2) {
  beta <- split(fit$models[[m]]$beta, rep(seq_along(fit$lnbeta[[m]]), fit$lnbeta[[m]]))
  predicted <- sapply(beta, function(b)
    as.vector(outer((age_grid - 24) / 12, 0:(length(b) - 1), `^`) %*% b))
  colors <- hcl.colors(ncol(predicted), "Dark 3")
  matplot(age_grid, predicted, type = "l", lty = 1, lwd = 2, col = colors,
          xlab = "Age (months)", ylab = indicators[m], xaxt = "n")
  axis(1, at = months)
  legend("topright", legend = paste0("G", seq_len(ncol(predicted))), col = colors, lwd = 2, bty = "n")
}
par(old_par)

