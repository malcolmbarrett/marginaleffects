#' @include get_coef.R
#' @rdname get_coef
#' @export
get_coef.ipw <- function(model, ...) {
    # An `ipw` result presents one of two readings, and which one is a property
    # of the object. `coef()` already returns the active one: the reported
    # causal effects in marginal mode, the weighted outcome model's
    # coefficients in conditional mode.
    stats::coef(model)
}


#' @include get_vcov.R
#' @rdname get_vcov
#' @export
get_vcov.ipw <- function(model, vcov = NULL, ...) {
    if (isFALSE(vcov)) {
        return(NULL)
    }
    if (isTRUE(checkmate::check_matrix(vcov))) {
        return(vcov)
    }

    # Both readings are covariances the fitting package computed over its own
    # stacked estimating equations. There is nothing here for `sandwich` or
    # `clubSandwich` to reweight, so a robust or clustered request is refused
    # rather than silently answered with the default. NULL and TRUE both mean
    # "the model's own covariance", and TRUE is what the main entry points pass.
    if (!is.null(vcov) && !isTRUE(vcov)) {
        stop_sprintf(
            "The `vcov` argument is not supported for models of class `ipw`. Its covariance is computed by the fitting package over the stacked estimating equations, and already accounts for the estimation of the weights. Supply a matrix to override it, or leave `vcov` at its default."
        )
    }

    stats::vcov(model)
}


#' @include set_coef.R
#' @rdname set_coef
#' @export
set_coef.ipw <- function(model, coefs, ...) {
    # Write back into whichever store the active reading presents, or the
    # jacobian comes back zero and every standard error is silently NA.
    #
    # The `[]` on the conditional branch is load bearing: it writes in place and
    # keeps the outcome model's coefficient names. Replacing the vector outright
    # drops them, and a perturbed model with unnamed coefficients loses its term
    # labels, which breaks any hypothesis string that names a term.
    if (identical(model[["effects"]], "conditional")) {
        model[["outcome_mod"]][["coefficients"]][] <- unname(coefs)
    } else {
        model[["estimates"]][["estimate"]] <- unname(coefs)
    }
    return(model)
}


#' @include get_predict.R
#' @rdname get_predict
#' @export
get_predict.ipw <- function(model, newdata, type = "response", ...) {
    # Predictions are a conditional-mode quantity. `sanitize_model_specific()`
    # refuses earlier with a better message, but a direct call must not read the
    # wrong thing either.
    if (!identical(model[["effects"]], "conditional")) {
        stop_sprintf(
            "Predictions require the conditional reading of an `ipw` result. Call `as_conditional()` on it first."
        )
    }
    pred <- stats::predict(
        model[["outcome_mod"]],
        newdata = newdata,
        type = type
    )
    rowid <- if (is.null(newdata[["rowid"]])) {
        seq_len(nrow(newdata))
    } else {
        newdata[["rowid"]]
    }
    out <- data.frame(rowid = rowid, estimate = unname(as.numeric(pred)))
    return(out)
}


#' @include sanity_model.R
#' @rdname sanitize_model_specific
#' @export
sanitize_model_specific.ipw <- function(
    model,
    calling_function = "unknown",
    ...
) {
    # `hypotheses()` is correct in either reading and on every standard error
    # method, so it is always let through. Every prediction-based entry point
    # arrives as "predictions" or "comparisons", slopes included.
    if (identical(calling_function, "hypotheses")) {
        return(model)
    }

    if (!identical(model[["effects"]], "conditional")) {
        stop_sprintf(
            "Predictions, comparisons, and slopes are conditional-mode quantities: they are per-row functions of the weighted outcome model, and this `ipw` result is presenting its marginal (population-averaged) effects. Call `as_conditional()` on it first, or use `hypotheses()` to test the marginal effects directly."
        )
    }

    # Some standard error methods have no stacked estimating system, so the
    # fitting package leaves the outcome model plain and its conditional reading
    # records no covariance. Prediction-based standard errors would then ignore
    # the estimation of the weights.
    if (!inherits(model[["outcome_mod"]], "ipw_model")) {
        stop_sprintf(
            "The conditional reading of this `ipw` result records no covariance from the joint estimation of the weights and the outcome, so standard errors here would ignore the estimation of the weights. Refit with a standard error method that produces one, or use `hypotheses()` on the marginal reading, which is correct either way."
        )
    }

    return(model)
}
