#' Cluster-based permutation tests for time series data, based on mixed-effects models or other \code{buildmer} models.
#' @param formula A normal formula, possibly using \code{lme4}-style random effects. This can also be a buildmer terms object, provided \code{dep} is passed in \code{buildmerControl}. Only a single response variable is supported. For binomial models, the \code{cbind} syntax is not supported; please convert your dependent variable to a proportion and use weights instead.
#' @param family The family.
#' @param data The data.
#' @template weightsoffset
#' @param series.var A one-sided formula giving the variable grouping the time series.
#' @template buildmer1
#' @param parallel Whether to parallelize the permutation testing using plyr's \code{parallel} option. Needs some additional set-up; see the plyr documentation.
#' @param progress A plyr \code{.progress} bar name, see the plyr documentation. If not \code{'none'} while \code{parallel=TRUE}, an ad-hoc solution will be used, which will be visible if the cluster nodes were created with \code{outfile=''}.
#' @template buildmer2
#' @examples
#' \donttest{
#' # Testing a single EEG electrode, with random effects by participants
#' perms <- clusterperm.lmer(Fz ~ Deviant * Session + (Deviant * Session | Subject),
#' 	data=MMN,series.var=~Time)
#' # Testing a single EEG electrode, with random effects by participants, ANOVA inference
#' perms <- clusterperm.lmer(Fz ~ Deviant * Session + (Deviant * Session | Subject),
#' 	data=MMN,series.var=~Time,type='anova')
#' }
#' \dontshow{
#' perms <- clusterperm.lmer(Fz ~ Deviant*Session + (1|Subject),data=MMN[MMN$Time > 200 & MMN$Time < 205,],series.var=~Time,nperm=2,type='anova')
#' perms <- clusterperm.lmer(Fz ~ Deviant*Session + (1|Subject),data=MMN[MMN$Time > 200 & MMN$Time < 205,],series.var=~Time,nperm=2,type='regression')
#' perms <- clusterperm.lmer(Fz ~ Session + (1|ppn),data=within(MMN[MMN$Time > 200 & MMN$Time < 205,],{Session <- factor(Session)}),series.var=~Time,nperm=2,type='regression')
#' }
#' @importFrom stats gaussian
#' @export
clusterperm.lmer <- function (formula,data=NULL,family=gaussian(),weights=NULL,offset=NULL,series.var=~0,buildmerControl=list(direction='order',crit='LRT',quiet=TRUE,ddf='lme4'),nperm=1000,type='regression',parallel=FALSE,progress='none') {
	return(clusterperm.work(buildmer::buildmer,formula,data,family,weights,offset,series.var,buildmerControl,nperm,type,parallel,progress))
}

#' Cluster-based permutation tests for time series data, based on \code{buildglmmTMB}.
#' @param formula A normal formula, possibly using \code{lme4}-style random effects. This can also be a buildmer terms object, provided \code{dep} is passed in \code{buildmerControl}. Only a single response variable is supported. For binomial models, the \code{cbind} syntax is not supported; please convert your dependent variable to a proportion and use weights instead.
#' @param family The family.
#' @param data The data.
#' @template weightsoffset
#' @param series.var A one-sided formula giving the variable grouping the time series.
#' @template buildmer1
#' @param parallel Whether to parallelize the permutation testing using plyr's \code{parallel} option. Needs some additional set-up; see the plyr documentation.
#' @param progress A plyr \code{.progress} bar name, see the plyr documentation. If not \code{'none'} while \code{parallel=TRUE}, an ad-hoc solution will be used, which will be visible if the cluster nodes were created with \code{outfile=''}.
#' @template buildmer2
#' @details
#' \code{clusterperm.glmmTMB} is much slower than \code{clusterperm.lmer}, but it is also more flexible, allowing for things like beta regression and zero-inflation.
#' @examples
#' \donttest{
#' # buildglmmTMB is much slower than clusterperm.lmer
#' \dontshow{if (FALSE)}
#' perms <- clusterperm.glmmTMB(Fz ~ Deviant * Session + (Deviant * Session | Subject),
#' 	data=MMN[MMN$Time > 150 & MMN$Time < 250,],series.var=~Time)
#' \dontshow{
#' perms <- clusterperm.glmmTMB(Fz ~ Deviant*Session + (1|Subject),data=MMN[MMN$Time > 200 & MMN$Time < 205,],series.var=~Time,nperm=2,type='regression')
#' }
#' }
#' @importFrom stats gaussian
#' @export
clusterperm.glmmTMB <- function (formula,data=NULL,family=gaussian(),weights=NULL,offset=NULL,series.var=~0,buildmerControl=list(direction='order',crit='LRT',quiet=TRUE),nperm=1000,type='regression',parallel=FALSE,progress='none') {
	if (type == 'anova') {
		stop('ANOVA is not available for glmmTMB models')
	}
	pkgcheck(c('glmmTMB','buildmer'))
	return(clusterperm.work(buildmer::buildglmmTMB,formula,data,family,weights,offset,series.var,buildmerControl,nperm,type,parallel,progress))
}

clusterperm.work <- function (buildmer,formula,data,family,weights,offset,series.var,buildmerControl,nperm,type,parallel,progress) {
	if (length(type) != 1 || !type %in% c('anova','regression')) {
		stop("Invalid 'type' argument (specify one of 'anova' or 'regression')")
	}
	dep <- if ('dep' %in% names(buildmerControl)) buildmerControl$dep else as.character(formula[2])
	if (all(is.null(weights))) {
		weights <- rep(1,length(data[[dep]]))
	}
	if (all(is.null(offset))) {
		offset <- rep(0,length(data[[dep]]))
	}
	if (is.null(data)) {
		y <- get(dep)
		x <- rep(0,length(y))
		data <- data.frame(.weights=x+1,.offset=x)
	} else {
		ix <- !is.finite(data[[dep]]) & !is.finite(weights) & !is.finite(offset)
		data <- data[ix,]
		if ('.weights' %in% names(data)) {
			stop("Please remove/rename the column named '.weights' from your data; this column name is used internally by permutes")
		}
		if ('.offset' %in% names(data)) {
			stop("Please remove/rename the column named '.offset' from your data; this column name is used internally by permutes")
		}
		data$.weights <- weights[ix]
		data$.offset <- offset[ix]
	}
	if (!inherits(series.var,'formula')) {
		series.var <- stats::reformulate(series.var)
	}
	series.var <- attr(terms(series.var),'term.labels')
	has.series <- length(series.var)
	if (has.series == 0) {
		timepoints <- rep(0,sum(ix))
	} else {
		if (has.series != 1) {
			stop('series.var does not appear to contain exactly one variable')
		}
		timepoints <- data[[series.var]]
		if (is.null(timepoints)) {
			stop('series.var ',series.var,' not found in data')
		}
	}
	check <- 'buildmer'
	if (has.series) {
		check <- c(check,'permuco')
	}
	if (type == 'anova') {
		check <- c(check,'car')
	}
	pkgcheck(check)

	if (is.character(family)) {
		family <- get(family)
	}
	if (is.function(family)) {
		family <- family()
	}
	if (has.series) {
		if (isTRUE(parallel)) {
			verbose <- progress != 'none'
			progress <- 'none'
		} else {
			if (!is.character(progress)) {
				stop("Invalid 'progress' specified for non-ad-hoc parallel solution (it has to be a character string)")
			}
			verbose <- FALSE
		}
	} else {
		verbose <- if (is.logical(progress)) progress else progress != 'none'
		progress <- 'none'
	}
	wrap <- function (t,fun,buildmer,formula,data,family,timepoints,buildmerControl,nperm,type,verbose) {
		errfun <- function (e) {
			# error in permutation-test function, return an empty result for this timepoint
			warning(e)
			data.frame(Factor=NA,p=NA)
		}
		ix <- timepoints == t
		data <- data[ix,]
		tryCatch(fun(t,buildmer,formula,data,family,timepoints,buildmerControl,nperm,type,verbose),error=errfun)
	}
	results <- plyr::alply(sort(unique(timepoints)),1,wrap,fit.buildmer,buildmer,formula,data,family,timepoints,buildmerControl,nperm,type,verbose,.parallel=parallel,.progress=progress,.inform=FALSE)
	terms <- lapply(results,`[[`,'terms')
	perms <- lapply(results,`[[`,'perms')
	df <- plyr::ldply(results,`[[`,'df',.id=series.var)
	for (i in seq_along(perms)) {
		# alply drops names for some very strange reason, so we saved them in the 'terms' element and restore them here
		names(perms[[i]]) <- terms[[i]]
	}

	df$p <- df$cluster_mass <- df$cluster <- rep(NA,nrow(df))
	# We need to invert the double-nested list from perm$time$factor to perm$factor$time
	for (x in unique(df$Factor)) {
		this.factor <- lapply(perms,`[[`,x) #all timepoints for this one factor
		df.LRT <- max(sapply(this.factor,function (x) x$df),na.rm=TRUE) #these will all be the same (because they are the same model comparison and these are ndf), except possibly in cases of rank-deficiency, hence why max is correct
		thresh <- stats::qchisq(.95,df.LRT) + .001 #+.001 to account for epsilon errors
		samp   <- sapply(this.factor,function (x) c(x$LRT,x$perms)) #columns are time, rows are samples
		p      <- apply(samp,2,function (x) if (sum(!is.na(x)) == 1) NA else mean(x[-1]+.001 >= x[1],na.rm=TRUE))
		# see GH issue #4: computing the cluster-mass test will fail in some cases, e.g. when only the intercept is involved (which is not permuted)
		stat   <- rep(NA,NCOL(samp)) #account for the possible failure case
		if (has.series) {
			cmass <- try(suppressWarnings(permuco::compute_clustermass(samp,thresh,sum,'greater'))$main,silent=TRUE)
			if (!inherits(cmass,'try-error')) {
				stat <- cmass #succeeded!
			}
		}
		df[df$Factor == x,c('p','cluster_mass','p.cluster_mass','cluster')] <- c(p,stat)
	}

	if (has.series && nrow(df)) {
		df <- cbind(df[,1],Measure=as.character(dep),df[,-1])
		colnames(df)[1] <- series.var
	} else {
		df[,1] <- df$cluster <- df$cluster_mass <- df$p.cluster_mass <- NULL
	}
	attr(df,'permutations') <- results
	class(df) <- c('permutes','data.frame')
	return(df)
}

fit.buildmer <- function (t,fun,formula,data,family,timepoints,buildmerControl,nperm,type,verbose) {
	buildmerControl$direction <- 'order'
	if (is.null(buildmerControl$quiet)) {
		buildmerControl$quiet <- TRUE
	}

	# First, make sure we only work with the tabular representation of the formula
	if (inherits(formula,'formula')) {
		if (is.null(buildmerControl$dep)) {
			buildmerControl$dep <- as.character(formula[2])
		}
		formula <- buildmer::tabulate.formula(formula)
	}
	fixed <- is.na(formula$grouping)
	terms <- stats::setNames(,formula$term[fixed])

	# For the fixed effects, we need to normalize the term/factor names (GH issue #4)
	# Don't try to be clever and set new.names <- old.names; we don't want to overwrite e.g. factors in the user data with model-matrix columns -- we don't know how they interact with the user's r.e. specification
	formula.fx <- buildmer::build.formula(NULL,formula[fixed,])
	mm <- stats::model.matrix(formula.fx,data)
	if (type == 'regression') {
		# one column per coefficient
		old.names <- colnames(mm)
		#terms <- apply(mm,2,identity,simplify=FALSE) #R >=4.1
		terms <- lapply(1:NCOL(mm),function (i) mm[,i])
	} else {
		# multiple columns per coefficient
		assign <- attr(mm,'assign')
		tl <- attr(terms(formula.fx),'term.labels')
		if (0 %in% assign) {
			tl <- c('(Intercept)',tl)
			assign <- assign + 1
		}
		old.names <- tl
		terms <- lapply(1:length(tl),function (i) {
			cols <- mm[,assign == i,drop=FALSE]
			colnames(cols) <- rep('',ncol(cols))
			cols
		})
	}
	new.names <- paste0('X',1:length(old.names))
	for (i in 1:length(new.names)) {
		if (new.names[i] %in% names(data)) {
			stop('Please remove column ',new.names[i],' from your data; its name conflicts with a name used internally in permutes!')
		}
		data[[new.names[i]]] <- terms[[i]]
	}
	new.formula <- data.frame(index=NA,grouping=NA,term=new.names,code=new.names,block=new.names)
	e <- environment(formula)
	formula <- rbind(new.formula,formula[!fixed,])
	environment(formula) <- e

	.weights <- data$.weights; .offset <- data$.offset #silence R CMD check warning
	buildmerControl$args <- c(buildmerControl$args,list(weights=.weights,offset=.offset))
	bm <- fun(formula=formula,data=data,family=family,buildmerControl=buildmerControl)
	# Apparently, on M1 systems, *something* can go wrong here (discovered by CRAN testing). In that case, return an empty model to avoid crashing in the logic aggregating the results
	if (inherits(try(model <- bm@model),'try-error')) {
		warning('Error fitting model:',bm)
		df <- data.frame(Factor=character(0),df=numeric(0),LRT=numeric(0),F=numeric(0))
		return(list(terms=character(0),perms=list(),df=numeric(0)))
	}

	perms <- lapply(1:length(new.names),function (i) {
		if (verbose) {
			time <- Sys.time()
			nmodels <- length(new.names) * length(unique(timepoints))
		}

		# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3883440/ propose that, to test a random effect:
		# 1. Get the marginal errors calculated by the alternative model
		# 2. To account for these errors' non-independence: weight the errors by the inverse of the random-effects correlation matrix
			# V0 = sigma^2_b_1_0 * ZtZ + sigma^2_e_0I
		# 3. Weight by Ut0^-1, where U0 = chol(V0)
			# lme4 parameterizes sigma^2_b_1_0 Z = Z Lambda_theta and we really want their transpose
		# 4. Permute the unweighted errors, then reweight the permuted data
		# 5. Reestimate *both* models with the fixed effects removed, which is necessary if any random effects happened to be similar

		# Our problem is much simpler than the above, however, since we only test fixed effects.
		# For us, the marginal errors are therefore based only on the fixed effects, and therefore *are* exchangeable under the null.
		# According to https://amstat.tandfonline.com/doi/abs/10.1080/01621459.1994.10476890, these marginal errors indeed are exchangeable even across subjects and items.

		# 1. Get the marginal errors based on quantities from the alternative model
		# These are y - XB - Zu, with the effect of interest *removed* from X
		# It's easier for us to just take the residuals and add this effect back in
		if (inherits(model,'glmmTMB')) {
			X <- lme4::getME(model,'X')
			B <- lme4::fixef(model)$cond
		} else if (inherits(model,'merMod')) {
			X <- lme4::getME(model,'X')
			B <- lme4::fixef(model)
		} else {
			X <- stats::model.matrix(formula(model),data)
			B <- stats::coef(model)
		}
		keep <- colnames(X) == new.names[i]
		X[,!keep] <- 0 #must modify X rather than B because X will be used as predictors in the permuted models
		if (any(wh <- !is.finite(B))) { #rank-deficiency in lm
			X[,wh] <- 0
		}
		e <- as.vector(stats::resid(model) + X %*% B)
		e <- family$linkinv(e)

		# 2/3. Random effects have already been partialed out, so these are independent and exchangeable
		# 4/5. Permute them and estimate a null and alternative model on the permuted data
		# The offset has been partialed out already, so will be ignored
		fit <- function (formula,data) {
			perm_weights <- data$.weights #silence R CMD check warning
			suppressWarnings(fun(formula,data,family=family,buildmerControl=list(REML=FALSE,quiet=TRUE,direction=NULL,calc.summary=FALSE,calc.anova=FALSE,args=list(weights=perm_weights))))@model
		}
		data <- list(
			perm_y = e,
			perm_X = X[,keep],
			perm_weights = data$.weights
		)
		# Optimization: nothing to do if the actually-kept columns are constant
		if (length(unique(as.vector(X[,keep]))) == 1) {
			perms <- NA
		} else {
			perms <- lapply(1:nperm,function (i) try({
				s <- sample(seq_along(e))
				data$perm_y <- data$perm_y[s]
				data$perm_weights <- data$perm_weights[s]
				m1 <- fit(perm_y ~ 0+perm_X,data)
				m0 <- fit(perm_y ~ 0,data)
				as.numeric(2*(stats::logLik(m1)-stats::logLik(m0)))
			},silent=TRUE))
			bad <- sapply(perms,inherits,'try-error')
			if (any(bad)) {
				perms[bad] <- NA
			}
			perms <- unlist(perms)
		}

		# Wrap up
		m1  <- fit(perm_y ~ 0+perm_X,data)
		m0  <- fit(perm_y ~ 0,data)
		ll1 <- stats::logLik(m1)
		ll0 <- stats::logLik(m0)
		LRT <- as.numeric(2*(ll1-ll0))
		df  <- attr(ll1,'df') - attr(ll0,'df')
		if (verbose) {
			diff <- Sys.time() - time
			cat('One of the',nmodels,'permutation tests finished in',diff,attr(diff,'units'),'\n')
		}
		list(perms=perms,LRT=LRT,df=df)
	})
	LRT <- sapply(perms,`[[`,'LRT')
	names(LRT) <- new.names

	scale.est <- !family(model)$family %in% c('binomial','poisson')
	is.mer <- inherits(model,'merMod')
	if (type == 'regression') {
		if (inherits(model,'glmmTMB')) {
			beta <- lme4::fixef(model)$cond
			vcov <- stats::vcov(model)$cond
		} else {
			beta <- (if (inherits(model,'merMod')) lme4::fixef else stats::coef)(model)
			vcov <- stats::vcov(model)
		}
		vcov  <- vcov[new.names,new.names]
		se    <- sqrt(diag(as.matrix(vcov))) #as.matrix needed to work around 'Error in diag(vcov(model)) : long vectors not supported yet: array.c:2186'
		LRT   <- LRT[new.names]
		beta  <- beta[new.names]
		se    <- se[new.names]
		df    <- data.frame(Factor=old.names,LRT=unname(LRT),beta=unname(beta),t=unname(beta/se))
		colnames(df)[4] <- if (scale.est) 't' else 'z'
		list(terms=old.names,perms=perms,df=df)
	} else {
		if (inherits(model,'gam')) {
			anovatab <- stats::anova(model) #is Type III
			Fvals <- anovatab$pTerms.chi.sq / anovatab$pTerms.df
			Fname <- 'F'
			df <- anovatab$pTerms.df
			names(df) <- names(Fvals) <- rownames(anovatab$pTerms.table)
		} else {
			if (is.mer) {
				anovatab <- car::Anova(model,type=3,test='Chisq')
				if (scale.est) {
					Fvals <- anovatab$Chisq / anovatab$Df
					Fname <- 'F'
				} else {
					Fvals <- anovatab$Chisq
					Fname <- 'Chisq'
				}
			} else {
				if (scale.est) {
					anovatab <- car::Anova(model,type=3,test='F')
					Fvals <- anovatab$F
					Fname <- 'F'
				} else {
					anovatab <- car::Anova(model,type=3,test='Wald')
					Fvals <- anovatab$Chisq
					Fname <- 'Chisq'
				}
				if (rownames(anovatab)[nr <- length(Fvals)] == 'Residuals') {
					anovatab <- anovatab[-nr,]
					Fvals <- Fvals[-nr]
				}
			}
			df <- anovatab$Df
			names(df) <- names(Fvals) <- rownames(anovatab)
		}
		df    <- df[new.names]
		LRT   <- LRT[new.names]
		Fvals <- Fvals[new.names]
		df <- data.frame(Factor=old.names,df=df,LRT=unname(LRT),F=unname(Fvals))
		colnames(df)[4] <- Fname
		list(terms=old.names,perms=perms,df=df)
	}
}
