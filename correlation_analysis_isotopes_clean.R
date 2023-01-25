# Import data (Note: the path needs to be set to point to the data)
setwd("/Users/imogenbrowne/Dropbox/Imogen/GRL paper/9_21_GRL/Stats Analyses/git_repository/")
isotopes_df = read.csv("BF_isotopes.csv")
plot(isotopes_df$Depth.mcd, isotopes_df$d18O.merged)
plot(isotopes_df$Depth.mcd, isotopes_df$d13C.merged)

# Compute Pearson's Correlation Coefficient
cor.test(isotopes_df$d18O.merged, isotopes_df$d13C.merged, 
         method = "pearson")


