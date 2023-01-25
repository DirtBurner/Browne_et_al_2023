# Import data (Note: the path needs to be set to point to the data)
# setwd("~/Projects/Side_Projects/imogen/wap")
setwd("/Users/imogenbrowne/Dropbox/Imogen/GRL paper/9_21_GRL/Stats Analyses/git_repository/")
TEX86_df = read.csv("TEX86_final.csv")
plot(TEX86_df$Depth.cmcd, TEX86_df$TEX86)
# plot(TEX86_df$Age.CE, TEX86_df$TEX86)

#Analytical uncertainty (pooled standard deviation of replicates) should replace those numbers that have no listed replicate uncertainty.
analytical_unc = 0.012 
TEX86_df$TEX86.stdev = replace(TEX86_df$TEX86.stdev, is.na(TEX86_df$TEX86.stdev), analytical_unc)

# Manually subsets the data to a given range
# curr_data_df = subset(TEX86_df,(TEX86_df$Age_CE>=1300) & 
#                         (TEX86_df$Age_CE<=1940))
curr_data_df = subset(TEX86_df,(TEX86_df$Depth.cmcd>19) &
                        (TEX86_df$Depth.cmcd<=100))
# curr_data_df = TEX86_df

# -----------------------------------------------------------------
# Plot data
plot(curr_data_df$Depth.cmcd, curr_data_df$TEX86, ylim = c(0.25,0.38))
arrows(curr_data_df$Depth.cmcd, curr_data_df$TEX86-curr_data_df$TEX86.stdev, 
       curr_data_df$Depth.cmcd, curr_data_df$TEX86+curr_data_df$TEX86.stdev, 
       length=0.05, angle=90, code=3)

# -----------------------------------------------------------------
# Weighted regression
#define weights to use
wt = 1/(curr_data_df$TEX86.stdev) # invert so that values with smaller uncertainty are weighted more.

#perform weighted least squares regression
wls_model = lm(TEX86 ~ Depth.cmcd, data = curr_data_df, weights = wt)

#view summary of model
summary(wls_model)

