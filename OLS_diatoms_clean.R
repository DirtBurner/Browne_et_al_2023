# Import data (Note: the path needs to be set to point to the data)
setwd("/Users/imogenbrowne/Dropbox/Imogen/GRL paper/9_21_GRL/Stats Analyses/git_repository/")
sea_ice_df = read.csv("diatom_sea_ice.csv")
plot(sea_ice_df$Age, sea_ice_df$F..curta.F..kerg)

# Manually subsets the data to a given range
curr_data_df = subset(sea_ice_df,(sea_ice_df$Age>1300) & 
                        (sea_ice_df$Age<1940))

# Plot current data
plot(curr_data_df$Age,curr_data_df$F..curta.F..kerg)

# Run the linear regression
curr_regression = lm(F..curta.F..kerg ~ Age, curr_data_df)
summary(curr_regression) # Look at results


