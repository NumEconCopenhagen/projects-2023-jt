
# Importing modules:
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import ipywidgets as widgets
from matplotlib_venn import venn2
import pandas_datareader.data as web
import plotly.express as px
from types import SimpleNamespace

# If you dont have the eurostat extension installed, it will be necessary to run the code below
## %pip install eurostat
import eurostat

class GDP_CapitaClass : 
    df = eurostat.get_data_df('nama_10_gdp')

    # We choose which rows that we want to see.
    # we have chosen to see the gross domestic product in Chain linked volumes (2015), million euro. 
    gdp = df[df['na_item'] == 'B1GQ']
    gdp = gdp[gdp['unit']=='CLV15_MEUR']
    
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

class PopulationClass :

    code = 'DEMO_PJAN'
    pars = eurostat.get_pars(code)
    df = eurostat.get_data_df(code, flags=False)

    # We access the data that we need
    # This time we are doing it by filtering directly from the data source instead of cleaning it afterwards.
    my_filter_pars = {'startPeriod':2012,'endPeriod': 2022, 'sex': 'T', 'age':'TOTAL'}
    population = eurostat.get_data_df(code, filter_pars=my_filter_pars)

    # We rename the column geo\TIME_PERIOD
    population.rename(columns={'geo\TIME_PERIOD': 'Country_code'}, inplace=True)

    # We are again deleting the columns we dont need. We are going by the aggregates and we therefore dont need these columns
    del_coloumns = ['freq' , 'unit', 'age', 'sex']
    population.drop(columns=del_coloumns, axis=1, inplace=True) 


class MergeClass :

    gdp = GDP_CapitaClass()
    population = PopulationClass()

    # We are now chaning the direction of the two datasets, making them long rather than wide. 
    population_long = pd.wide_to_long(population , stubnames='' , i= 'Country_code', j= 'year')

    gdp_long = pd.wide_to_long(gdp, stubnames= '', i= 'Country_code' , j= 'year')

    # We will now merge the two datasets, by doing an inner join; meaning we choose the observations (countries) which are in both datasets. 
    inner = pd.merge(gdp_long, population_long, how = 'inner' , on = ['Country_code' , 'year'])

    # We are now renaming the columns
    inner.rename(columns={'_x':'GDP', '_y':'Population'}, inplace=True)

    # Dropping the countries that have Nan for all values of either GDP or population
    inner.dropna(inplace=True)
        
    # We are now resetting the index
    inner.reset_index(inplace = True)

    # We are now creating a new column for GDP per capita, since GDP is in millions we multiply by 1.000.000
    inner["GDP_Cap"] = inner["GDP"]*1000000/inner["Population"]

    # We are now merging the data with the excel file
    df = pd.read_excel('C_Name_ISO3.xlsx')

    # We are now merging the two datasets
    merge = pd.merge(inner, df, how = 'inner', on = ['Country_code'])

    # we will now reaarange the columns
    merge = merge[[ 'Country_Name','Country_code', 'ISO_3_Code', 'year', 'GDP', 'Population', 'GDP_Cap']]
        

class Plot_class :
    def __init__(self):
        par = self.par = SimpleNamespace()
        par.merge = MergeClass.merge()       
    
    def plot_choropleth(self) : 
        par = self.par
        merge = par.merge

        # We will delete kosovo as the ISO code is not available to use for the choropleth map
        merge = merge[merge['ISO_3_Code'] != 'XKX']
        # We are now creating the choropleth map
        fig = px.choropleth(merge, 
                            locations='ISO_3_Code',
                            color = 'GDP_Cap',
                            scope = 'europe',
                            hover_name='Country_Name',
                            animation_frame='year')

        return fig.show()    
      
    def line_interactive(self) :
        par = self.par
        merge = par.merge
    
        def plot_line(Country_Name) :
            #We start by creating the normal line plot

            I = merge['Country_Name'] == Country_Name
            # We are now setting the x-axis o show the years and the y-axis to show the GDP/Cap
            ax = merge.loc[I,:].plot(x='year', y='GDP_Cap', style='-o', legend=False) 

            # We are now setting the title and the labels
            ax.set_title('GDP/Capita 2012-2022 for '+ Country_Name)
            ax.set_ylabel('GDP/Capita in euros')
            ax.set_xlabel('Year')
            ax.set_xlim(merge['year'].min(), merge['year'].max())
            ax.set_xticks(np.arange(merge['year'].min(), merge['year'].max()+1))
        
            return plt.show() 
        # Now we make the lineplot interactive

        line = widgets.interact(plot_line, 
                            inner = widgets.fixed(merge), 
                            Country_Name = widgets.Dropdown( name = 'Country_Name',
                                                            options = merge['Country_Name'].unique(),
                                                            value = 'Denmark'))

        return line        
     
    
    def scatter_interactive(self) :
    
        par = self.par
        merge = par.merge

        def plot_scatter(year) : 
            # We are now creating the scatterplot

            I = merge['year'] == year
            ax = merge.loc[I,:].plot(x='GDP_Cap', y='Population', style='o', legend=False)
            ax.set_ylabel('Population in millions')
            ax.set_xlabel('GDP per capita in euros')
            ax.set_title(f"Scatterplot of GDP per capita and Population for {year}")
            plt.subplots_adjust(left=0.2, right=1, top=0.9, bottom=0.1)
    
            return plt.show 

        # Now we make the scatterplot interactive
        year_widget = widgets.Dropdown(options=merge['year'].unique(), value=2022, description='Year:')
        scatter = widgets.interact(plot_scatter, inner=widgets.fixed(merge), year=year_widget)

        return scatter
        
    