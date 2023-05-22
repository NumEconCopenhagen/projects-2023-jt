
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

class GDP_CapitaClass : 

    def __init__(self):
        pass

    def Get_GDP(self):

        df = eurostat.get_data_df('nama_10_gdp')

        # We choose which rows that we want to see.
            # we have chosen to see the gross domestic product in Chain linked volumes (2015), million euro. 
        gdp = df[df['na_item'] == 'B1GQ']
        gdp = gdp[gdp['unit']=='CLV15_MEUR']

        return gdp
    
    def Clean_GDP (self) : 
            
        gdp = self.Get_GDP() 

        # We remove the columns freq, unit, na_item, and the years 1975-2011
        drop_these = ['freq' ,] + [str(i) for i in range(1975,2012,1)]
        gdp.drop(drop_these, axis=1, inplace=True)

        # We rename the coloumn geo\TIME_PERIOD
        gdp.rename(columns={'geo\TIME_PERIOD': 'Country_code'}, inplace=True)

        # We remove the aggregate values in our data as we only are interested in the specific countries.
        remove_these = ['EA', 'EA12', 'EA19', 'EA20', 'EU15', 'EU27_2020', 'EU28']

        for i in remove_these : 
            gdp = gdp[gdp['Country_code']!= i]

        # we are resetting the index
        gdp.reset_index(inplace = True, drop = True)

        return gdp
    

    def Get_Population(self, do_print=False):
        
        code = 'DEMO_PJAN'
        pars = eurostat.get_pars(code)
        df = eurostat.get_data_df(code, flags=False)

        # We access the data that we need
        # This time we are doing it by filtering directly from the data source instead of cleaning it afterwards.
        my_filter_pars = {'startPeriod':2012,'endPeriod': 2022, 'sex': 'T', 'age':'TOTAL'}
        population = eurostat.get_data_df(code, filter_pars=my_filter_pars)

        return population    
    
    def Clean_Population(self, do_print=False) :
         
        population = self.Get_Population()

        # We rename the column geo\TIME_PERIOD
        population.rename(columns={'geo\TIME_PERIOD': 'Country_code'}, inplace=True)

        # We are again deleting the columns we dont need. We are going by the aggregates and we therefore dont need these columns
        del_coloumns = ['freq' , 'unit', 'age', 'sex']

        population.drop(columns=del_coloumns, axis=1, inplace=True) 

        return population
    
    def Merge_Data(self) :
        gdp = self.Clean_GDP()
        population = self.Clean_Population()

        # We are now chaning the direction of the two datasets, making them long rather than wide. 
        population_long = pd.wide_to_long(population , stubnames='' , i= 'Country_code', j= 'year')

        gdp_long = pd.wide_to_long(gdp, stubnames= '', i= 'Country_code' , j= 'year')

        # We will now merge the two datasets, by doing an inner join; meaning we choose the observations (countries) which are in both datasets. 
        inner = pd.merge(gdp_long, population_long, how = 'inner' , on = ['Country_code' , 'year'])

        return inner

    def Clean_merge(self) : 
        inner = self.Merge_Data()

        # We are now renaming the columns
        inner.rename(columns={'_x':'GDP', '_y':'Population'}, inplace=True)

        # Dropping the countries that have Nan for all values of either GDP or population
        inner.dropna(inplace=True)
        
        # We are now resetting the index
        inner.reset_index(inplace = True)

        # We are now creating a new column for GDP per capita, since GDP is in millions we multiply by 1.000.000
        inner["GDP/Cap"] = inner["GDP"]*1000000/inner["Population"]

        return inner
    
    def Merge_excel(self):
        inner = self.Clean_merge()

        # We are now merging the data with the excel file
        df = pd.read_excel('C_Name_ISO3.xlsx')

        # We are now merging the two datasets
        merge = pd.merge(inner, df, how = 'inner', on = ['Country_code'])

        # we will now reaarange the columns
        merge = merge[[ 'Country_name','Country_code', 'ISO3', 'year', 'GDP', 'Population', 'GDP/Cap']]
        return merge