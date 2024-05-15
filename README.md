Code and data from Walsworth et al. (in press at _Ecology of Freshwater Fish_) "Interactions between runoff volume, timing, and annual temperatures shape migration phenology of a threatened adfluvial sucker"

File descriptions:
- DOYConversion.csv: CSV file with data to convert Julian dates to calendar dates - used in figures only
- JS_MigrationTiming_EnvironmentalData.csv: CSV file containing environmental data and principal components used in linear regression analyses.
- JS_MigrationTiming_EnvironmnetalRegressions.R: R script running linear regressions of environmental predictors against annual MCMC posterior estimates of peak migration date, among individual variation, and residence time.
- JS_MigrationTiming_Model_Gamma_JAGScode.bug: JAGS code to run the June sucker migration timing model assuming individual migration days are gamma distributed (not included in paper due to poor relative fit, included here for optional comparison).
- JS_MigrationTiming_Model_Normal_JAGScode.bug: JAGS code to run the June sucker migration timing model assuming individual migration days are normally distributed (the focal model in the paper).
- JS_MigrationTiming_Model_Posteriors.rdata: RData file containing JAGS output for normal migration model. Included MCMC posterior distributions of all parameters. Needed for regression analyses.
- JS_MigrationTiming_RunJAGSModel.R: R script to read in data and run the migration timing model through JAGS.
