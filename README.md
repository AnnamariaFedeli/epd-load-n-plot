![](https://i.imgur.com/VbfxACB.png)

## epd_plot.py

This Python script is for plotting and extracting data from the raw data of Solar Orbiters EPT, HET and STEP instruments. <br>
It works in conjunction with epd_loader, which fetches and loads the raw data from [ESA](http://soar.esac.esa.int/soar).

### Requirements

* [epd_loader.py](https://github.com/jgieseler/solo-loader), download and place in same folder

### Usage

Plotting, for the most part, follows a three-step formula:
1. Load raw data via epd_loader.
2. Run data through extract_data -function.
3. Use returned data in plot functions.

Following is a usage example for loading and plotting EPT data:

```python
from epd_plot import *

start_date = 20201210
end_date = 20201211

# Load ept data.
df_protons, df_electrons, energies = read_epd_cdf('ept', 'sun', 'l2', start_date, end_date, path='~/solo/data/', autodownload=True)

plot_start = '2020-12-10-2100'
plot_end = '2020-12-11-0200'
search_start = '2020-12-10-2330'
search_end = '2020-12-11-0010'
bg_start = '2020-12-10-2200'
bg_end = '2020-12-10-2300'
averaging_mode = 'mean'
averaging = 10

# extract_data function returns an object that can be used in plotting functions as input.
data = extract_data(df_protons, df_electrons, plot_start, plot_end, search_start, search_end, bg_start, bg_end, instrument='ept', data_type='l2', averaging_mode=averaging_mode, averaging=averaging)

# Plots all available energy channels in a single image.
plot_channels(data, bg_subtraction=False)

# Plots the spectrum of energy channels.
plot_spectrum(data, bg_subtraction=True)
```

Plotting HET data works in a similar way:

```python
from epd_plot import *

start_date = 20201210
end_date = 20201211

# Load het data.
df_protons, df_electrons, energies = read_epd_cdf('het', 'sun', 'l2', start_date, end_date, path='~/solo/data/', autodownload=True)

plot_start = '2020-12-10-2100'
plot_end = '2020-12-11-0200'
search_start = '2020-12-10-2330'
search_end = '2020-12-11-0010'
bg_start = '2020-12-10-2200'
bg_end = '2020-12-10-2300'
averaging_mode = 'mean'
averaging = 10

# extract_data function returns an object that can be used in plotting functions as input.
data = extract_data(df_protons, df_electrons, plot_start, plot_end, search_start, search_end, bg_start, bg_end, instrument='het', data_type='l2', averaging_mode=averaging_mode, averaging=averaging)

# Plots all available energy channels in a single image.
plot_channels(data, bg_subtraction=False)

# Plots the spectrum of energy channels.
plot_spectrum(data, bg_subtraction=True)
```
Plotting STEP data works much in the same way, with a few differences: epd_loader only returns a single dataframe, that is used in extract_step_data -function, instead of the extract_data -function. There's also a masking -parameter for STEP data in extract_step_data, which toggles ion contamination masking.

```python
from epd_plot import *

start_date = 20201210
end_date = 20201211

# Load step data.
df_particles, energies = epd_load(sensor='step', viewing='sun', level='l2', startdate=start_date, enddate=end_date, path='~/solo/data/', autodownload=True)

plot_start = '2020-12-10-2000'
plot_end = '2020-12-11-0200'
search_start = '2020-12-10-2330'
search_end = '2020-12-11-0030'
bg_start = '2020-12-10-2030'
bg_end = '2020-12-10-2230'
averaging_mode = 'mean'
averaging = 10

# extract_data function returns an object that can be used in plotting functions as input.
data = extract_step_data(df_particles, plot_start, plot_end, search_start, search_end, bg_start, bg_end, instrument='step', data_type='l2', averaging_mode=averaging_mode, averaging=averaging, masking=True)

# Plots all available energy channels in a single image.
plot_channels(data, bg_subtraction=False)

# Plots the spectrum of energy channels.
plot_spectrum(data, bg_subtraction=True)
```

An example Jupyter Notebook is provided, which can be used to quickly make and save plots and data.

### Functions

**extract_data**(df_protons, df_electrons, plotstart, plotend, searchstart, searchend, bgstart, bgend, *instrument*, *data_type*, *averaging_mode*, *averaging*, *masking*, *ion_conta_corr*)

This function extracts the data required to use the other functions from provided proton and electron data. This is the sort of 'core' function of the program that is always ran first.

Arguments:
- df_protons: The proton dataframe provided by epd_loader.py.
- df_electrons: The electron dataframe provided by epd_loader.py.
- plotstart: Start time of the plot. Format: yyyy-mm-dd-hhmm
- plotend: End time of the plot
- searchstart: Start time of the peak search period.
- searchend: End time of the peak search period.
- bgstart: Start time of the background period.
- bgend: End time of the background period.
- *instrument*: (default: 'ept') (options: 'ept', 'het', 'step') Provide the instrument in question.
- *data_type*: (default: 'L2') (options: 'l2', 'll') Provide level of data.
- *averaging_mode*: (default: 'none') (options: 'none', 'mean', 'rolling_window') Averaging mode. 'rolling_window' is not recommended.
- *averaging*: (default: 2) Provide averaging time in minutes, if 'mean' toggled. For 'rolling_window' this is the window size in observations.
- *masking*: (default: False) Toggled the STEP ion masking.
- *ion_conta_corr*: (default: False) Toggled the ion contamination correction for EPT.
- *frac_nan_threashold*: (default 0.4) Used to to check if there is enough non-nan flux data points in the search-period interval. 
        If not, the flux and uncertainty value of that energy channel are set to nan and therefore excluded from the spectrum.
---
**plot_channels**(data, *bg_subtraction*, *savefig*, *path*, *key*)

Plots the electron flux channels of provided data.

Arguments:
- data: This is the return provided by extract_data or extract_step_data function.
- *bg_subtraction*: (default: False) Applies background subtraction to the channel plot.
- *savefig*: (default: False) Toggles saving of image file.
- *path*: (default: '') Provide a path for saved image, if toggled.
- *key*: (default: '') This string is added to the end of the saved image filename.

---
**plot_spectrum**(data, *bg_subtraction*, *savefig*, *path*, *key*)

Plots the spectrum plot of provided data.

Arguments:
- data: This is the return provided by extract_data or extract_step_data function.
- *bg_subtraction*: (default: True) Applies background subtraction to the spectrum plot.
- *savefig*: (default: False) Toggles saving of image file.
- *path*: (default: '') Provide a path for saved image, if toggled.
- *key*: (default: '') This string is added to the end of the saved image filename.

---
**plot_check**(data)

Unfinished function that plots all flux channels in the same figure and colors them.

Arguments:
- data: The return provided by either extract_data or extract_step_data.

---
**write_to_csv**(data, *path*, *key*)

Function to write relevant data to a human-readable .csv -format.

Arguments:
- data: The return provided by either extract_data or extract_step_data.
- *path*: (default: '') Provide a path for the saved .csv -file.
- *key*: (default: '') This string is added to the end of the saved .csv -file filename. 

---
**acc_flux**(data, *time*)

Unfinished function for calculating the accumulated electron flux from a given time period.

Arguments:
- data: The return provided by either extract_data or extract_step_data.
- *time*: (default: []) Provide a time period from which to compute the flux. Format: [yyyy-mm-dd-hhmm, yyyy-mm-dd-hhmm]. This default to the search period if left empty.
