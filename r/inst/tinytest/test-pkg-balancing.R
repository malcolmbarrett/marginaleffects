source("helpers.R")
using("marginaleffects")
if (!requiet("balancing")) {
    exit_file("balancing not installed")
}
if (!requiet("causalgenerics")) {
    exit_file("causalgenerics not installed")
}

set.seed(1024)
N <- 400
x1 <- rnorm(N)
exposure <- rbinom(N, 1, plogis(0.8 * x1))
y <- rbinom(N, 1, plogis(-0.5 + 0.8 * exposure + 0.9 * x1))
dat <- data.frame(x1 = x1, exposure = exposure, y = y)

# balancing solves for the weights directly rather than modelling a propensity
# score, and builds the same `ipw` class from them. `exposure_type` is passed
# explicitly because balance() reports the type it infers otherwise.
bw <- balancing::balance(
    dat,
    exposure,
    x1,
    method = balancing::bw_entropy(),
    estimand = "ate",
    exposure_type = "binary"
)
dat$w <- stats::weights(bw)
# quasibinomial() solves the same estimating equation as binomial() and does not
# warn that non-integer weights are not counts.
mod <- balancing::ipw(
    bw,
    glm(y ~ exposure + x1, data = dat, family = quasibinomial(), weights = w)
)

# The default reading is marginal: the reported causal effects. The labels come
# from the fit itself rather than from a hardcoded list.
known <- as.data.frame(mod)
expect_equal(mod$effects, "marginal")
expect_equal(names(coef(mod)), known$term)

# hypotheses() tests whichever reading the object presents. as.data.frame() on
# an `ipw` result reports the tidier column names: term, estimate, std.error,
# statistic, p.value.
h <- hypotheses(mod)
expect_equal(h$term, known$term)
# `expect_equivalent` wherever the two sides need not agree on attributes and
# names: marginaleffects estimate columns may carry a label, and the named
# vectors further down have names the other side does not.
expect_equivalent(h$estimate, known$estimate, tolerance = 1e-10)
expect_equal(h$std.error, known$std.error, tolerance = 1e-6)

# Per-row quantities are refused in the marginal reading, with the remedy named.
expect_error(
    avg_comparisons(mod, variables = "exposure", newdata = dat),
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

# Contrasts and slopes on the conditional reading are the same calls made
# against the wrapped outcome model directly, which needs no support for the
# `ipw` class.
cmp <- avg_comparisons(cond, variables = "exposure", newdata = dat)
cmp_route_a <- avg_comparisons(mod$outcome_mod, variables = "exposure", newdata = dat)
expect_equivalent(cmp$estimate, cmp_route_a$estimate, tolerance = 1e-10)
expect_equal(cmp$std.error, cmp_route_a$std.error, tolerance = 1e-6)

slo <- avg_slopes(cond, variables = "x1", newdata = dat)
slo_route_a <- avg_slopes(mod$outcome_mod, variables = "x1", newdata = dat)
expect_equivalent(slo$estimate, slo_route_a$estimate, tolerance = 1e-10)
expect_equal(slo$std.error, slo_route_a$std.error, tolerance = 1e-6)

# Averaged contrasts on the conditional reading recover the reported risk
# difference, because an equally weighted average of unit-level contrasts is the
# ATE this fit targets.
expect_equivalent(
    cmp$estimate,
    known$estimate[known$term == "rd"],
    tolerance = 1e-6
)

# hypotheses() on the conditional reading tests outcome terms against the
# corrected covariance.
hc <- hypotheses(cond)
expect_equal(hc$term, names(coef(mod$outcome_mod)))
expect_equivalent(hc$std.error, sqrt(diag(vcov(cond))), tolerance = 1e-6)

# Joint tests read one parameter vector with its own covariance, so they are
# correct without needing a guard.
j <- hypotheses(cond, joint = TRUE)
expect_true(is.finite(j$statistic))

# A fit targeting a focal estimand standardizes to a different population, so
# the standardization weights have to be supplied. balancing has no propensity
# score to tilt by: its weights are design quantities, and the counterpart for a
# focal estimand is the indicator for the focal group.
bw_att <- balancing::balance(
    dat,
    exposure,
    x1,
    method = balancing::bw_entropy(),
    estimand = "att",
    exposure_type = "binary"
)
dat$w_att <- stats::weights(bw_att)
att_fit <- balancing::ipw(
    bw_att,
    glm(y ~ exposure + x1, data = dat, family = quasibinomial(), weights = w_att)
)
known_att <- as.data.frame(att_fit)
att_cond <- causalgenerics::as_conditional(att_fit)
att_focal <- avg_comparisons(
    att_cond,
    variables = "exposure",
    newdata = dat,
    wts = as.numeric(dat$exposure == 1)
)
expect_equivalent(
    att_focal$estimate,
    known_att$estimate[known_att$term == "rd"],
    tolerance = 1e-6
)
# Omitting the weights averages over the wrong population. The size of the gap
# is a property of the fixture, so only its presence is asserted.
att_flat <- avg_comparisons(att_cond, variables = "exposure", newdata = dat)
expect_true(
    abs(att_flat$estimate - known_att$estimate[known_att$term == "rd"]) > 1e-6
)

# A robust variance request has nothing to act on and is refused.
expect_error(
    avg_comparisons(cond, variables = "exposure", newdata = dat, vcov = "HC3"),
    pattern = "not supported for models of class"
)
