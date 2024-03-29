pkgcheck <- function (needPkgs) {
	missing <- !sapply(needPkgs,requireNamespace)
	if (nmissing <- sum(missing)) {
		if (nmissing == 1) {
			stop(paste('Please install the following missing package:',needPkgs[missing]))
		} else {
			stop(paste('Please install the following missing packages:',needPkgs[missing],collapse=', '))
		}
	}
	if ('buildmer' %in% needPkgs && utils::packageVersion('buildmer') < '2.3') {
		stop('Package buildmer must be at least version 2.3; please update')
	}
}
