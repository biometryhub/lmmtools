##' Generalised Heritability for a single trial
##'
##' @rdname genherit1
##' @export
##'
genherit1.asreml <- function(asr, id = 'Genotype', only = NULL, Gmat = NULL, ...) {
    stopifnot(requireNamespace("asreml"))
    response <- as.character(asr$call[[2]][2])
    if(is.null(only)) {
        if(!is.null(Gmat)) stop("\n Error:  argument only is NULL, but Gmat is not NULL\n")
        only <- id
    }
    if(packageVersion("asreml") == "3.0") {
        mypred <- asreml::predict.asreml(asr, classify=id, maxiter=1, only=only, vcov=TRUE, ...)
        which.vc <- grep(id, names(asr$gammas))
    } else {
        mypred <- asreml::predict.asreml(asr, classify=id, maxit=1, only=only, vcov=TRUE, ...)
        which.vc <- grep(id, names(asr$vparameters), fixed=TRUE)
    }
    my.vcov <- mypred$vcov
    if(is.null(Gmat)) {
##################  version 3 might require summary(asr)$varcomp[which.vc,2]
        h2 <- 1- sum(diag(my.vcov))/(asreml::summary.asreml(asr)$varcomp[which.vc,1]*nrow(my.vcov))
    }  else {
        if(is.null(attr(Gmat, "INVERSE"))) stop('\n Gmat must have an "INVERSE" attribute\n')
        if(!attr(Gmat, "INVERSE")) Gmat <- MASS::ginv(Gmat)
        h2 <- 1- sum(diag(Gmat%*%my.vcov))/nrow(my.vcov)
    }
    names(h2) <- response
    h2
}
