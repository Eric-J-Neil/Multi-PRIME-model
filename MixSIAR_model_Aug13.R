# Plant Water Uptake Program


library(MixSIAR)

# Load mix data
mix.filename <- system.file("MixSIAR_Aug13_Consumer.csv", package = "MixSIAR")
mix <- load_mix_data(filename=mix.filename,
					 iso_names=c("H","O"),
					 factors=NULL,
					 fac_random=NULL,
					 fac_nested=NULL,
					 cont_effects=NULL)

# Load source data
source.filename <- system.file("MixSIAR_Aug13_Source.csv", package = "MixSIAR")
source <- load_source_data(filename=source.filename,
						   source_factors=NULL,
						   conc_dep=FALSE,
						   data_type="raw",
						   mix)

# Load discrimination/TDF data
discr.filename <- system.file("MixSIAR_Aug13_Discrimination.csv", package = "MixSIAR")
discr <- load_discr_data(filename=discr.filename, mix)


#############################################################################################
# # INFORMATIVE prior (construct alpha from Volumetric Water Content and Root Length Density)
#############################################################################################


# My Prior Information
kw.alpha <- c(0.04,2.06,0.98,0.92)

# Generate alpha hyperparameters scaling sum(alpha)=n.sources
kw.alpha <- kw.alpha*length(kw.alpha)/sum(kw.alpha)

# the Dirichlet hyperparameters for the alpha.prior cannot be 0 (but can set = .01)
kw.alpha[which(kw.alpha==0)] <- 0.01

setwd("C:/Users/Admin/Documents/R/win-library/4.1/...")

# Make isospace plot
plot_data(filename="isospace_plot",
		  plot_save_pdf=FALSE,
		  plot_save_png=TRUE,
		  mix,source,discr)

dev.off()

# Plot your informative prior
plot_prior(alpha.prior=kw.alpha,
		   source=source,
		   plot_save_pdf=TRUE,
		   plot_save_png=FALSE,
		   filename="prior_plot_kw_inf")

dev.off()

###################################### Run Model - very long

setwd("C:/Users/Admin/Documents/R/win-library/4.1/...")

# Define model structure and write JAGS model file
model_filename <- "MixSIAR_model_kw_inf.txt"   # Name of the JAGS model file
resid_err <- TRUE
process_err <- FALSE
write_JAGS_model(model_filename, resid_err, mix, source)

 jags.inf <- run_model(run="very long",mix,source,discr,model_filename,alpha.prior=kw.alpha)
 output_JAGS(jags.inf, mix, source)

dev.off()
dev.off()
dev.off()

