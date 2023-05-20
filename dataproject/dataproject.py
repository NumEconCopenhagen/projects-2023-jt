def keep_regs(df, regs):
    """ Example function. Keep only the subset regs of regions in data.

    Args:
        df (pd.DataFrame): pandas dataframe 

    Returns:
        df (pd.DataFrame): pandas dataframe

    """ 
    
    for r in regs:
        I = df.reg.str.contains(r)
        df = df.loc[I == False] # keep everything else
    
    return df
# Importing modules:
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import ipywidgets as widgets
from matplotlib_venn import venn2
import pandas_datareader.data as web

# If you dont have the eurostat extension installed, it will be necessary to run the code below
## %pip install eurostat
import eurostat

# We access the data from eurostat and call it df
df = eurostat.get_data_df('nama_10_gdp')

# We choose which rows that we want to see.
    # we have chosen to se the gross domestic product in Chain linked volumes (2015), million euro. 
gdp = df[df['na_item'] == 'B1GQ']
gdp = gdp[gdp['unit']=='CLV15_MEUR']

# We remove the columns freq, unit, na_item, and the years 1975-2011
drop_these = ['freq' ,] + [str(i) for i in range(1975,2012,1)]
gdp.drop(drop_these, axis=1, inplace=True) # axis = 1 -> columns, inplace=True -> changed, no copy made

 # We rename the coloumn geo\TIME_PERIOD
gdp.rename(columns={'geo\TIME_PERIOD': 'Country_code'}, inplace=True)

# We remove the aggregate values in our data as we only are interested in the specific countries.
remove_these = ['EA', 'EA12', 'EA19', 'EA20', 'EU15', 'EU27_2020', 'EU28']

for i in remove_these : 
    gdp = gdp[gdp['Country_code']!= i]

# we are resetting the index
gdp.reset_index(inplace = True, drop = True)




# We are now downloading our data for the population of the different countries.
# we name our parameters
code = 'DEMO_PJAN'
pars = eurostat.get_pars(code)

# We access the data that we need
# This time we are doing it by filtering directly from the data source instead of cleaning it afterwards.
my_filter_pars = {'startPeriod':2012,'endPeriod': 2022, 'sex': 'T', 'age':'TOTAL'}
population = eurostat.get_data_df(code, filter_pars=my_filter_pars)

# We rename the column geo\TIME_PERIOD

population.rename(columns={'geo\TIME_PERIOD': 'Country_code'}, inplace=True)

# We are again deleting the columns we dont need. We are going by the aggregates and we therefore dont need these columns
del_coloumns = ['freq' , 'unit', 'age', 'sex']

population.drop(columns=del_coloumns, axis=1, inplace=True) 

# We are now chaning the direction of the two datasets, making them long rather than wide. 
population_long = pd.wide_to_long(population , stubnames='' , i= 'Country_code', j= 'year')

gdp_long = pd.wide_to_long(gdp, stubnames= '', i= 'Country_code' , j= 'year')

# We will now merge the two datasets, by doing an inner join; meaning we choose the observations (countries) which are in both datasets. 
inner = pd.merge(gdp_long, population_long, how = 'inner' , on = ['Country_code' , 'year'])



# We are now renaming our columns to the representative names
inner.rename(columns={'_x': 'GDP', '_y':'Population'}, inplace=True)
inner.reset_index(inplace=True)

# Dropping the countries that have Nan for all values of either GDP or population
inner.dropna(inplace=True)

# We are now creating a new column for GDP per capita, since GDP is in millions we multiply by 1.000.000
inner["GDP/Cap"] = inner["GDP"]*1000000/inner["Population"]