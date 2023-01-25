# Import data (Note: the path needs to be set to point to the data)
setwd("/Users/imogenbrowne/Dropbox/Imogen/GRL paper/9_21_GRL/Stats Analyses/git_repository/")

forams_df = read.csv("BF_isotopes.csv")
plot(forams_df$Depth.mcd, forams_df$d13C.merged)
plot(forams_df$Depth.mcd,forams_df$d18O.merged)

#Analytical uncertainty (pooled standard deviation of replicates) should replace those numbers that have no listed replicate uncertainty.

analytical_unc_d18O = 0.094 
forams_df$X18O.stdev = replace(forams_df$X18O.stdev, is.na(forams_df$X18O.stdev), analytical_unc_d18O)
analytical_unc_d13C = 0.153
forams_df$d13C.stdev = replace(forams_df$d13C.stdev, is.na(forams_df$d13C.stdev), analytical_unc_d13C)

# Manually subsets the data to a given range
curr_data_df = subset(forams_df,(forams_df$Depth.mcd>0) &
                        (forams_df$Depth.mcd<=.44))

# Plot data d18O
plot(curr_data_df$Depth.mcd, curr_data_df$d18O.merged, ylim = c(3,3.8))
arrows(curr_data_df$Depth.mcd, curr_data_df$d18O.merged-curr_data_df$X18O.stdev, 
       curr_data_df$Depth.mcd, curr_data_df$d18O.merged+curr_data_df$X18O.stdev, 
       length=0.05, angle=90, code=3)

#Plot data d13C
plot(curr_data_df$Depth.mcd, curr_data_df$d13C.merged, ylim = c(-0.8,0.5))
arrows(curr_data_df$Depth.mcd, curr_data_df$d13C.merged-curr_data_df$d13C.stdev, 
       curr_data_df$Depth.mcd, curr_data_df$d13C.merged+curr_data_df$d13C.stdev, 
       length=0.05, angle=90, code=3)


# Run the linear regression (OLS; non-weighted)
curr_regression_d13C = lm(d13C.merged ~ Depth.mcd, curr_data_df)
summary(curr_regression_d13C)

curr_regression_d18O = lm(d18O.merged ~ Depth.mcd, curr_data_df)
summary(curr_regression_d18O)

# -----------------------------------------------------------------
# Weighted regression just to check what they look like
#define weights to use
wt_d13C = 1/(curr_data_df$d13C.stdev) # invert so that values with smaller uncertainty are weighted more.

wt_d18O = 1/(curr_data_df$X18O.stdev)

#perform weighted least squares regression
wls_model_d13C = lm(d13C.merged ~ Depth.mcd, data = curr_data_df, weights = wt_d13C)
wls_model_d18O = lm(d18O.merged ~ Depth.mcd, data = curr_data_df, weights = wt_d18O)

#view summary of model
summary(wls_model_d13C)
summary(wls_model_d18O)