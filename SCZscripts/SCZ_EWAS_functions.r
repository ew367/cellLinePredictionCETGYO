

# Sex EWAS
Sex <- function(row, pheno){
	full.model <- lm(row ~ Sex + Age + Phenotype + Tissue.Centre, data=QCmetrics)
	stats <- coef(summary(full.model))['SexM',c('Estimate','Std. Error','Pr(>|t|)')]
	return(stats)
}

# Sex*Phenotype EWAS
SexGroup <- function(row, pheno){
	full.model <- lm(row ~ Sex*Phenotype + Age + Sex + Phenotype + Tissue.Centre, data=QCmetrics)
	stats <- c(coef(summary(full.model))['SexM:PhenotypeSchizophrenia',c('Estimate','Std. Error','Pr(>|t|)')], coef(summary(full.model))['SexM',c('Estimate','Std. Error','Pr(>|t|)')], coef(summary(full.model))['PhenotypeSchizophrenia',c('Estimate','Std. Error','Pr(>|t|)')])
	return(stats)
}
