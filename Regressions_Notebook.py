# %% [markdown]
# # TEX86 Regressions - Antarctic Peninsula
#
# This workbook goes through regression techniques to answer the question of whether the trend in time of reconstructed TEX86 temperatures is significant or not. The question was posed by the lone reviewer after the initial submission of this manuscript. To answer it, we factor the analytical and *geological* uncertainty of each TEX86 measurement to weight the data for regression. We then employ two different modes of regression. 
#

# %%
#Load packages and data, set x, y, and uncertainty variables, and print first 5 lines of data:
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

data_df = pd.read_csv('TEX_data.csv')
#Analytical uncertainty should replace those numbers that have no listed replicate uncertainty.
#According to the manuscript, that value is: 0.015 or use mean standard deviation of replicates.
analytical_unc = np.mean(data_df['TEX86_stdev'])
data_df['TEX86_stdev'].fillna(analytical_unc, inplace=True)
depth_mask = [19, 100]
data_df = data_df.mask(data_df.mask(data_df['Depth_cmcd']>=min(depth_mask))['Depth_cmcd']<=max(depth_mask)).dropna()

x = data_df['Depth_cmcd']
y = data_df['TEX86']
unc_y = data_df['TEX86_stdev']

print(data_df.head())

# %%
#Plot data with error bars

fig1, axes1 = plt.subplots(nrows=1, ncols=1)

axes1.errorbar(
    x,
    y,
    unc_y, 
    marker='o',
    markersize=10,
    mec='k',
    color='pink',
    ecolor='k',
)

#Add trend:
a, b = np.polyfit(x, y, deg=1)
model = a*x + b

axes1.plot(x, model, linestyle='--', color='k')
axes1.legend(['Data', 'OLS Trend (no weights)'])

plt.ylabel('TEX 86 Index')
plt.xlabel('Year CE')
print(a, b)


# %% [markdown]
# ## Assessing the Trend
#
# It is apparent that there is a trend, but how much do the individual uncertainties contribute to the trend? And how significant is the trend? We'll deal with the first question first.

# %%
#How much do the individual uncertaintes contribute to the trend?
print(a, b)
X = np.c_[np.ones(x.shape[0]),x]
BA = np.linalg.inv(X.T @ X) @ X.T @ y
print(BA)
#Linear algebraic method matches the single line of code that permitted the model in the figure above.
#Now what about weights?
#print(data_df['TEX86_stdev'])
unc = np.diag([1/val for val in unc_y])
BA_w = np.linalg.inv(X.T @ np.linalg.inv(unc) @ X) @ (X.T @ np.linalg.inv(unc) @ y)
print(BA_w)

model_w = BA_w[0] + BA_w[1]*x



# %%
#Plot data as above, but with new weighted model:

fig2, axes2 = plt.subplots(nrows=1, ncols=1)

axes2.errorbar(
    x,
    y,
    unc_y, 
    marker='o',
    markersize=10,
    mec='k',
    color='pink',
    ecolor='k',
)

#Add old trend:
model = a*x + b
axes2.plot(x, model, linestyle='--', color='grey')
#Add new trent:
axes2.plot(x, model_w, linestyle='dashdot', color='k')
axes2.legend(['OLS Trend (no weights)', 'OLS Trend (weighted)','Data' ])

plt.ylabel('TEX 86 Index')
plt.xlabel('Year CE')

plt.savefig('Browne_et_al_2023_OLS_Models.svg')

print(BA_w[0])


# %% [markdown]
# ## Shallower Trend, But Significant?
# The dash-dot trend above is through the points weighted by uncertainty. It is flattened (closer to 0) compared to the dashed trend. That means that, given the same statistics, the weighted trend is less likely to be significant (i.e. more likely to include a slope of 0 in the range of likely slopes). Below, we find out if a Monte Carlo approach can elucidate the probability of the slope containing 0 given thousands of permutations of the data within established bounds. 

# %%
#Development of functions to do Monte Carlo analysis

def regress_data(x, y):
    slope, intercept = np.polyfit(x, y, deg=1)
    return slope, intercept

def jiggle_data(y, unc_y):
    y_jig = [np.random.uniform(val-unc, val+unc) for val, unc in zip(y, unc_y)]
    return y_jig

def multi_OLS(x, y, unc_y, iterations):
    slope_list = []
    intercept_list = []
    model_list = []
    for v in range(iterations):
        s, i = regress_data(x, jiggle_data(y, unc_y))
        slope_list.append(s)
        intercept_list.append(i)
        mod = x*s+i
        model_list.append(mod)

    return(slope_list, intercept_list, model_list)



    
    




# %%
slopes, intercepts, models = multi_OLS(x, y, 2*unc_y, 1000)

fig3, axes3 = plt.subplots(nrows=1, ncols=1)
axes3.errorbar(
    x,
    y,
    unc_y, 
    marker='o',
    markersize=10,
    mec='k',
    color='pink',
    ecolor='k',
)
for ind, val in enumerate(slopes):
    axes3.plot(
        x,
        models[ind],
        color='grey',
        marker=None
    )

#Add old trend:
model = a*x + b
axes3.plot(x, model, linestyle='--', color='k')
#Add new trent:
axes3.plot(x, model_w, linestyle='dashdot', color='k')

plt.ylabel('TEX 86 Index')
plt.xlabel('Year CE')

plt.savefig('Browne_et_al_2023_AllModels.svg')


# %% [markdown]
# ## Comparison of Monte Carlo and Weighted Regressions
# Interestingly, for this data set, the weighted regression deviates from the pool of Monte Carlo models formed from "jiggling" the data points within 2x their standard deviation. Because these standard deviations encompass both analytical and geological uncertainty, the "jiggling" of the data serves to imitate a scientist resampling the core X times and re-measuring the values. 96% of those repeated measurements will lie within +/- 2 standard deviations of the value measured herein. 
#
# **Why does the error-weighted model deviate?** Well, the functionality of weighting for errors nullifies those data points with large uncertainties and the model stablizes with more trust in those without large uncertainties. On the other hand, the Monte Carlo approach treats all data points equally, but allows them to vary within uncertainties. Because there are quite a few data points that do not have large uncertainties, the model is quite a bit more stable than indicated by the difference between the error-weighted and unweighted regressions. 
#
# **But what about the question at hand? *Is the trend statistically significant?*** To answer this, let's take a look at a simple histogram:

# %%
#Plot histograms of the slopes and intercepts of the Monte Carlo approach

fig4, axes4 = plt.subplots(nrows=2, ncols=1)
axes4[0].hist(slopes, color='peru', edgecolor='k')
ylims=axes4[0].get_ylim()
axes4[0].vlines(a, min(ylims), max(ylims), color='k', linestyle='--')
axes4[0].vlines(BA[1], min(ylims), max(ylims), color='k', linestyle='dashdot')
axes4[0].set_title('Slopes')
axes4[1].hist(intercepts, color='mediumslateblue', edgecolor='k')
ylims=axes4[1].get_ylim()
axes4[1].vlines(b, min(ylims), max(ylims), color='k', linestyle='--')
axes4[1].vlines(BA[0], min(ylims), max(ylims), color='k', linestyle='dashdot')
axes4[1].set_title('Intercepts')

plt.savefig('Browne_et_al_2023_histograms.svg')
print('Average slope = ', np.mean(slopes), ' +/- ', np.std(slopes))
