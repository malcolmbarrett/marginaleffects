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

# `pool_ipw()` computes both readings on every call, presents the one the call
# named, and stores the other whole. The reading a pooled result presents is
# therefore a view of one pooling rather than something settled at pool time:
# `as_conditional()` and `as_marginal()` move a result between the two by
# exchanging which one it presents, and the accessors take a per-call `effects`
# argument that reports the other reading without changing the result. The two
# objects below are the same pooling presented two ways. Everything up to the
# flip expectations at the end of this file reads whichever reading its object
# presents; those expectations check the moves between them.
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

# Both readings come out of one pooling, so moving a pooled result to the other
# reading has to report what pooling directly into that reading reports. The
# flip exchanges stored frames and recomputes nothing, so the two agree to the
# last bit here; the tolerances are the ones this file uses elsewhere rather
# than a claim about how far apart they could fall.
hfc <- hypotheses(causalgenerics::as_conditional(pooled))
expect_equal(hfc$term, hc$term)
expect_equivalent(hfc$estimate, hc$estimate, tolerance = 1e-10)
expect_equal(hfc$std.error, hc$std.error, tolerance = 1e-6)

# The move in the other direction, which is a separate method rather than the
# one above read backwards.
hfm <- hypotheses(causalgenerics::as_marginal(pooled_cond))
expect_equal(hfm$term, h$term)
expect_equivalent(hfm$estimate, h$estimate, tolerance = 1e-10)
expect_equal(hfm$std.error, h$std.error, tolerance = 1e-6)

# Out to the other reading and back is the result that went in, since the flip
# exchanges the presented frames with the stored ones and exchanging twice
# restores both. Nothing is recomputed on the way, so the round trip is compared
# exactly rather than to a tolerance.
hrt <- hypotheses(
    causalgenerics::as_marginal(causalgenerics::as_conditional(pooled))
)
expect_identical(hrt$term, h$term)
expect_identical(as.numeric(hrt$estimate), as.numeric(h$estimate))
expect_identical(hrt$std.error, h$std.error)

# The `set_coef()` round trip has to survive the flip. `coef()` reads the
# presented frame's `estimate` column and `set_coef()` writes it, and a flip
# moves which frame that is, so a non-linear hypothesis on a flipped result is
# the expectation that would catch a perturbation written back to the reading no
# longer presented: the jacobian would come back zero and the standard error
# silently NA. The conditional reading's terms are the outcome models'
# coefficient names, so the transform is of the exposure coefficient `z`.
pooled_flipped <- causalgenerics::as_conditional(pooled)
known_flipped <- as.data.frame(pooled_flipped)
b_z <- known_flipped$estimate[known_flipped$term == "z"]
se_z <- known_flipped$std.error[known_flipped$term == "z"]
hz <- hypotheses(pooled_flipped, "exp(`z`) = 1")
expect_equivalent(hz$estimate, exp(b_z) - 1, tolerance = 1e-8)
expect_equal(hz$std.error, exp(b_z) * se_z, tolerance = 1e-6)

# The refusal is about what a pooled result holds rather than about which
# reading it presents, so flipping one does not turn per-row quantities back on.
expect_error(
    predictions(causalgenerics::as_conditional(pooled), newdata = head(dat, 3)),
    pattern = "retains no fitted outcome model"
)

# Naming a reading on an accessor reports that reading for the one call and
# leaves the result presenting the one it presented before, which is what every
# expectation above was written against.
cc <- coef(pooled, effects = "conditional")
expect_identical(names(cc), names(coef(pooled_cond)))
expect_equal(pooled$effects, "marginal")
expect_identical(names(coef(pooled)), known$term)
