source("helpers.R")
using("marginaleffects")
if (!requiet("propensity")) {
    exit_file("propensity not installed")
}
if (!requiet("causalgenerics")) {
    exit_file("causalgenerics not installed")
}

set.seed(1024)
N <- 400
x1 <- rnorm(N)
z <- rbinom(N, 1, plogis(0.8 * x1))
y <- rbinom(N, 1, plogis(-0.5 + 0.8 * z + 0.9 * x1))
dat <- data.frame(x1 = x1, z = z, y = y)

ps_mod <- glm(z ~ x1, data = dat, family = binomial())
dat$w <- propensity::wt_ate(ps_mod, .exposure = dat$z, exposure_type = "binary")
# quasibinomial() solves the same estimating equation as binomial() and does not
# warn that non-integer weights are not counts.
out_mod <- glm(y ~ z + x1, data = dat, family = quasibinomial(), weights = w)
mod <- propensity::ipw(ps_mod, out_mod, estimand = "ate")

# The default reading is marginal: the reported causal effects.
expect_equal(mod$effects, "marginal")
expect_equal(names(coef(mod)), c("rd", "log(rr)", "log(or)"))

# hypotheses() tests whichever reading the object presents. as.data.frame() on
# an `ipw` result reports the tidier column names: term, estimate, std.error,
# statistic, p.value.
known <- as.data.frame(mod)
h <- hypotheses(mod)
expect_equal(h$term, known$term)
# `expect_equivalent` where a marginaleffects column is one side of the
# comparison: `estimate` carries a `label` attribute, and the named vectors
# further down have names the other side does not.
expect_equivalent(h$estimate, known$estimate, tolerance = 1e-10)
expect_equal(h$std.error, known$std.error, tolerance = 1e-6)

# A non-linear transform, checked against the delta method by hand. A string
# hypothesis reports lhs - rhs, hence the "- 1".
rr <- exp(known$estimate[known$term == "log(rr)"])
se_lnrr <- known$std.error[known$term == "log(rr)"]
ht <- hypotheses(mod, "exp(`log(rr)`) = 1")
expect_equivalent(ht$estimate, rr - 1, tolerance = 1e-8)
expect_equal(ht$std.error, rr * se_lnrr, tolerance = 1e-6)

# Per-row quantities are refused in the marginal reading, with the remedy named.
expect_error(
    avg_comparisons(mod, variables = "z", newdata = dat),
    pattern = "conditional-mode quantities"
)
expect_error(
    predictions(mod, newdata = head(dat)),
    pattern = "conditional-mode quantities"
)

# The conditional reading presents the weighted outcome model.
cond <- causalgenerics::as_conditional(mod)
expect_equal(names(coef(cond)), names(coef(mod$outcome_mod)))
expect_equal(unname(coef(cond)), unname(coef(mod$outcome_mod)))

p <- predictions(cond, newdata = head(dat))
expect_equal(nrow(p), 6)
expect_equal(
    p$estimate,
    unname(predict(mod$outcome_mod, newdata = head(dat), type = "response")),
    tolerance = 1e-10
)

# Averaged contrasts on the conditional reading recover the reported risk
# difference, because an equally weighted average of unit-level contrasts is the
# ATE this fit targets.
cmp <- avg_comparisons(cond, variables = "z", newdata = dat)
expect_equal(
    cmp$estimate,
    known$estimate[known$term == "rd"],
    tolerance = 1e-6
)

# newdata is optional: causalgenerics registers model.frame.ipw, so
# insight::get_data() recovers the estimation data from the outcome model.
cmp_default <- avg_comparisons(cond, variables = "z")
expect_equal(cmp_default$estimate, cmp$estimate, tolerance = 0)
expect_equal(cmp_default$std.error, cmp$std.error, tolerance = 0)
# The weights column the outcome model's frame carries is excluded, or it would
# be read as a variable the model was fitted on. Pinned as the whole column set
# rather than as the absence of one name, so that a method dropping more than it
# should is caught too.
expect_equal(
    colnames(model.frame(mod)),
    setdiff(colnames(model.frame(mod$outcome_mod)), "(weights)")
)

# Slopes on the conditional reading are the same call made against the wrapped
# outcome model directly, which needs no support for the `ipw` class.
slo <- avg_slopes(cond, variables = "x1", newdata = dat)
slo_route_a <- avg_slopes(mod$outcome_mod, variables = "x1", newdata = dat)
expect_equivalent(slo$estimate, slo_route_a$estimate, tolerance = 1e-10)
expect_equal(slo$std.error, slo_route_a$std.error, tolerance = 1e-6)

# hypotheses() on the conditional reading tests outcome terms against the
# corrected covariance.
hc <- hypotheses(cond)
expect_equal(hc$term, names(coef(mod$outcome_mod)))
expect_equivalent(hc$std.error, sqrt(diag(vcov(cond))), tolerance = 1e-6)

# Joint tests read one parameter vector with its own covariance, so they are
# correct without needing a guard.
j <- hypotheses(cond, joint = TRUE)
expect_true(is.finite(j$statistic))

# Categorical exposure: one row per non-reference level.
if (!requiet("nnet")) {
    exit_file("nnet not installed")
}
set.seed(2048)
a <- factor(sample(c("a", "b", "c"), N, replace = TRUE))
yc <- rbinom(N, 1, plogis(-0.5 + 0.4 * (a == "b") + 0.8 * (a == "c") + 0.3 * x1))
dat_cat <- data.frame(x1 = x1, a = a, y = yc)
ps_cat <- nnet::multinom(a ~ x1, data = dat_cat, trace = FALSE)
dat_cat$w <- propensity::wt_ate(
    predict(ps_cat, type = "probs"),
    dat_cat$a,
    exposure_type = "categorical"
)
out_cat <- glm(y ~ a, data = dat_cat, family = quasibinomial(), weights = w)
mod_cat <- propensity::ipw(ps_cat, out_cat, estimand = "ate")

cmp_cat <- avg_comparisons(
    causalgenerics::as_conditional(mod_cat),
    variables = "a",
    newdata = dat_cat
)
expect_equal(nrow(cmp_cat), 2)
known_cat <- as.data.frame(mod_cat)
expect_equal(
    cmp_cat$estimate,
    known_cat$estimate[known_cat$term == "rd"],
    tolerance = 1e-4
)

# A fit whose outcome model records no corrected covariance is refused for the
# prediction surfaces, and allowed for hypotheses().
mod_lin <- propensity::ipw(
    ps_mod,
    glm(y ~ z, data = dat, family = quasibinomial(), weights = w),
    estimand = "ate",
    se_method = "linearization"
)
expect_error(
    avg_comparisons(
        causalgenerics::as_conditional(mod_lin),
        variables = "z",
        newdata = dat
    ),
    pattern = "records no covariance"
)
h_lin <- hypotheses(mod_lin)
expect_equal(h_lin$std.error, as.data.frame(mod_lin)$std.error, tolerance = 1e-6)

# A robust variance request has nothing to act on and is refused.
expect_error(
    avg_comparisons(cond, variables = "z", newdata = dat, vcov = "HC3"),
    pattern = "not supported for models of class"
)
