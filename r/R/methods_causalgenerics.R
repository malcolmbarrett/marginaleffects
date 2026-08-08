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


#' @rdname get_coef
#' @export
get_coef.ipw_pooled <- function(model, ...) {
    # An `ipw_pooled` result presents one of two readings, and which one was
    # settled when it was pooled. `coef()` already returns the active one: the
    # pooled causal effects in marginal mode, the pooled coefficients of the
    # weighted outcome models in conditional mode.
    stats::coef(model)
}


#' @rdname get_vcov
#' @export
get_vcov.ipw_pooled <- function(model, vcov = NULL, ...) {
    if (isFALSE(vcov)) {
        return(NULL)
    }
    if (isTRUE(checkmate::check_matrix(vcov))) {
        return(vcov)
    }

    # Both readings are pooled covariances: Rubin's rules combining covariances
    # the fitting package computed over its own stacked estimating equations.
    # There is nothing here for `sandwich` or `clubSandwich` to reweight, so a
    # robust or clustered request is refused rather than silently answered with
    # the default. NULL and TRUE both mean "the model's own covariance", and
    # TRUE is what the main entry points pass.
    if (!is.null(vcov) && !isTRUE(vcov)) {
        stop_sprintf(
            "The `vcov` argument is not supported for models of class `ipw_pooled`. Its covariance combines, by Rubin's rules, covariances the fitting package computed over the stacked estimating equations of each pooled analysis, so it already accounts for the estimation of the weights and for the variation between analyses. Supply a matrix to override it, or leave `vcov` at its default."
        )
    }

    stats::vcov(model)
}


#' @rdname set_coef
#' @export
set_coef.ipw_pooled <- function(model, coefs, ...) {
    # Write back into the store the active reading presents, or the jacobian
    # comes back zero and every standard error is silently NA. The presented
    # reading always lives in the `estimates` frame, and a flip swaps the
    # frames, so writing that frame's `estimate` column writes exactly where
    # `coef()` reads in either reading. The alternate reading stashed alongside
    # it is left stale on purpose: this perturbed copy is internal to the
    # delta-method pipeline, and is never flipped or read through a per-call
    # `effects` argument.
    model[["estimates"]][["estimate"]] <- unname(coefs)
    return(model)
}


#' @rdname sanitize_model_specific
#' @export
sanitize_model_specific.ipw_pooled <- function(
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

    # The same refusal in both readings, because the reason is the same in both
    # and flipping the reading is no remedy: `as_conditional()` changes which
    # reading a pooled result presents, not what it holds, and a pooled result
    # is pooled coefficients rather than a fit, so neither reading has per-row
    # quantities.
    stop_sprintf(
        "A pooled result reports pooled effects and retains no fitted outcome model, so per-row quantities such as predictions, comparisons, and slopes are unavailable in either reading. `hypotheses()` is the supported surface for an `ipw_pooled` result. To work with per-row quantities, call marginaleffects on each unpooled fit's conditional reading before pooling."
    )
}
