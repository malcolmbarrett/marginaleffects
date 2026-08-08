source("helpers.R")
using("marginaleffects")
if (!requiet("propensity")) {
    exit_file("propensity not installed")
}
if (!requiet("causalgenerics")) {
    exit_file("causalgenerics not installed")
}

# One analysis: a seeded dataset, the propensity model fitted to it, and the
# weighted outcome model those two make an `ipw` result from. Three of them
# stand in for three imputations of one incomplete dataset, which is what
# `pool_ipw()` pools. Independently simulated rather than perturbed copies of
# one dataset, so the between-imputation variance the pooling adds is not
# degenerate.
analyze_one <- function(seed, n = 400) {
    set.seed(seed)
    x1 <- rnorm(n)
    z <- rbinom(n, 1, plogis(0.8 * x1))
    y <- rbinom(n, 1, plogis(-0.5 + 0.8 * z + 0.9 * x1))
    d <- data.frame(x1 = x1, z = z, y = y)
    ps_mod <- glm(z ~ x1, data = d, family = binomial())
    d$w <- propensity::wt_ate(ps_mod, .exposure = d$z, exposure_type = "binary")
    # quasibinomial() solves the same estimating equation as binomial() and does
    # not warn that non-integer weights are not counts.
    out_mod <- glm(y ~ z + x1, data = d, family = quasibinomial(), weights = w)
    list(fit = propensity::ipw(ps_mod, out_mod, estimand = "ate"), data = d)
}

analyses <- lapply(c(1024, 2048, 4096), analyze_one)
fits <- lapply(analyses, function(x) x$fit)
dat <- analyses[[1]]$data

# `pool_ipw()` chooses the reading at pool time, and a pooled result carries no
# way to change it afterwards: there is no `as_conditional()` for this class, so
# the two readings are two objects rather than two views of one.
pooled <- causalgenerics::pool_ipw(fits)
pooled_cond <- causalgenerics::pool_ipw(fits, effects = "conditional")
expect_equal(pooled$effects, "marginal")
expect_equal(pooled_cond$effects, "conditional")

# The labels come from the pooled result itself rather than from a hardcoded
# list, and the two accessors have to agree on them, since `hypotheses()` reads
# the names off `coef()` and every expectation below names a row by its `term`.
known <- as.data.frame(pooled)
expect_identical(names(coef(pooled)), known$term)

# hypotheses() tests the reading the pooled result presents. as.data.frame() on
# an `ipw_pooled` result reports term, estimate, std.error, statistic, df, and
# p.value.
h <- hypotheses(pooled)
expect_equal(h$term, known$term)
# `expect_equivalent` where a marginaleffects column is one side of the
# comparison: `estimate` carries a `label` attribute the other side does not.
expect_equivalent(h$estimate, known$estimate, tolerance = 1e-10)
expect_equal(h$std.error, known$std.error, tolerance = 1e-6)
expect_equivalent(h$statistic, h$estimate / h$std.error, tolerance = 1e-10)

# The pooled result refers its statistic and its p-value to t on the row's own
# Barnard-Rubin degrees of freedom, and `hypotheses()` refers them to the normal
# unless `df` says otherwise, so the p-values above are deliberately not
# compared. `df` takes one number for the whole table rather than one per row: a
# vector is accepted only at the length of `newdata`, and a call made against a
# model has no rows of it. Every pooled effect has its own count, so the t-based
# inference is reproduced one effect at a time, by selecting a single row and
# naming that row's count.
for (i in seq_len(nrow(known))) {
    hd <- hypotheses(
        pooled,
        hypothesis = sprintf("`%s` = 0", known$term[i]),
        df = known$df[i]
    )
    expect_equal(nrow(hd), 1L)
    expect_equivalent(hd$p.value, known$p.value[i], tolerance = 1e-6)
}

# A non-linear transform of the pooled marginal reading, checked against the
# delta method by hand. A string hypothesis reports lhs - rhs, hence the "- 1".
# This is the expectation that exercises the `set_coef()` round trip: the
# perturbed estimates have to be written back into the store `coef()` reads
# from, or the jacobian comes back zero and the standard error is silently NA.
rr <- exp(known$estimate[known$term == "log(rr)"])
se_lnrr <- known$std.error[known$term == "log(rr)"]
ht <- hypotheses(pooled, "exp(`log(rr)`) = 1")
expect_equivalent(ht$estimate, rr - 1, tolerance = 1e-8)
expect_equal(ht$std.error, rr * se_lnrr, tolerance = 1e-6)

# The three marginal effects are three functions of the same two weighted
# marginal means, so each pooled analysis contributed a covariance that is
# singular in exact arithmetic, and a joint test over all three inverts a matrix
# that is non-singular here only because three imputations happened to
# contribute three different null spaces. Two effects are two functions of two
# means, which is the test that is well posed however the fixture falls.
j <- hypotheses(pooled, joint = c("rd", "log(rr)"))
expect_true(is.finite(j$statistic))

# The conditional reading pools the outcome models' coefficients, which are a
# full rank parameter vector, so the whole set is testable jointly.
jc <- hypotheses(pooled_cond, joint = TRUE)
expect_true(is.finite(jc$statistic))

# The conditional reading tests those coefficients against the pooled
# covariance, whose diagonal is the square of the standard errors the frame
# reports, since both are the same combination of the same numbers.
hc <- hypotheses(pooled_cond)
expect_equal(hc$term, names(coef(pooled_cond)))
expect_equivalent(hc$std.error, sqrt(diag(vcov(pooled_cond))), tolerance = 1e-6)

# A pooled result holds no component models at all: `pool_ipw()` keeps the
# pooled estimates, their covariance, and the pooling diagnostics, and nothing
# else. There is no outcome model to predict from, so per-row quantities are
# impossible by construction rather than merely unsupported, and both readings
# are refused for that one reason. Pooled coefficients are not a fitted model.
expect_error(
    predictions(pooled, newdata = head(dat, 3)),
    pattern = "retains no fitted outcome model"
)
expect_error(
    avg_comparisons(pooled, variables = "z", newdata = head(dat, 3)),
    pattern = "retains no fitted outcome model"
)
expect_error(
    predictions(pooled_cond, newdata = head(dat, 3)),
    pattern = "retains no fitted outcome model"
)
expect_error(
    avg_comparisons(pooled_cond, variables = "z", newdata = head(dat, 3)),
    pattern = "retains no fitted outcome model"
)

# A robust variance request has nothing to act on, exactly as on an unpooled
# result: the pooled covariance is Rubin's combination of covariances the
# fitting package computed over its own stacked estimating equations.
expect_error(
    hypotheses(pooled, vcov = "HC3"),
    pattern = "not supported for models of class"
)
