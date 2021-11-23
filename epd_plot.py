from epd_loader import *
from matplotlib.ticker import LinearLocator, MultipleLocator, AutoMinorLocator
import matplotlib.dates as mdates
from matplotlib.ticker import FormatStrFormatter
from adjustText import adjust_text

def extract_data(df_protons, df_electrons, plotstart, plotend, searchstart, searchend, bgstart, bgend, instrument = 'ept', data_type = 'l2', averaging_mode='none', averaging=2, masking=False, ion_conta_corr=False):
    """ determines an energy spectrum from time series data for any of the Solar Orbiter / EPD sensors

    Parameters
    ----------
    df_protons : pandas DataFrame
        contains proton (ion) data if instrument is 'het' ('ept')
    df_electrons : pandas DataFrame
        electron data
    plotstart : string
        start time of the time series plot, e.g., '2020-11-18-0000'
    plotend : string
        end time of the time series plot, e.g., '2020-11-18-2230'
    searchstart : string
        start time of the time interval used for determining the spectrum
    searchend : string
        end time of the time interval used for determining the spectrum
    bgstart : string
        start time of the background time interval (used for background subtraction)
    bgend : string
        end time of the background time interval (used for background subtraction)
    instrument : str, optional
        'ept', 'het', or 'step'; by default 'ept'
    data_type : str, optional
        which data level (e.g., low latency (ll) or level2 (l2)) is used. This affects the number of energy channels; by default 'l2'
    averaging_mode : str, optional
        averaging of the data, 'mean', 'rolling_window', or 'none'; by default 'none'
    averaging : int, optional
        number of minutes used for averaging, by default 2
    masking : bool, optional
        Refers only to STEP data. If true, time intervals with significant (5 sigma) ion contamination are masked; by default False
    ion_conta_corr : bool, optional
        Refers only to EPT data. If true, ion contamination correction is applied; by default False

    Returns
    -------
    df_electron_fluxes : pandas DataFrame
        data frame that cotains the electron flux data with applied averaging and optional contamination correction (EPT) or masking (STEP)
    df_info : pandas DataFrame
        data frame that contains the spectrum data and all its metadata (which is saved to csv in the function write_to_csv())
    [searchstart, searchend]: list of strings
    [e_low, e_high] : list of float
    [instrument, data_type] . list of strings
    """
    # Takes proton and electron flux and uncertainty values from original data.
    if(instrument != 'step'):

        df_electron_fluxes = df_electrons['Electron_Flux'][plotstart:plotend]
        df_electron_uncertainties = df_electrons['Electron_Uncertainty'][plotstart:plotend]

    if(instrument == 'ept'):

        df_proton_fluxes = df_protons['Ion_Flux'][plotstart:plotend]
        df_proton_uncertainties = df_protons['Ion_Uncertainty'][plotstart:plotend]

        if(data_type == 'll'):

            channels = [0,1,2,3,4,5,6,7]

            for i in channels:
                df_electron_fluxes = df_electron_fluxes.rename(columns={'Ele_Flux_{}'.format(i):'Electron_Flux_{}'.format(i)})
                df_electron_uncertainties = df_electron_uncertainties.rename(columns={'Ele_Flux_Sigma_{}'.format(i):'Electron_Uncertainty_{}'.format(i)})

            e_low = [0.0329, 0.0411, 0.0537, 0.0733, 0.1013, 0.1425, 0.1997, 0.2821]
            e_high = [0.0411, 0.0537, 0.0733, 0.1013, 0.1425, 0.1997, 0.2821, 0.3977]

        elif(data_type == 'l2'):

            channels = range(0,34) 

            e_low = [0.0312, 0.0330, 0.0348, 0.0380, 0.0406, 0.0432, 0.0459, 0.0497, 0.0533, 0.0580, 0.0627, 0.0673, 0.0731, 0.0788, 0.0856, 0.0934, 0.1011, 0.1109, 0.1197, 0.1305, 0.1423, 0.1541, 0.1679, 0.1835, 0.1995, 0.2181, 0.2371, 0.2578, 0.2817, 0.3061, 0.3339, 0.3661, 0.3989, 0.4348]
            e_high = [0.0348, 0.0369, 0.0380, 0.0406, 0.0432, 0.0459, 0.0497, 0.0533, 0.0580, 0.0627, 0.0673, 0.0731, 0.0788, 0.0856, 0.0934, 0.1011, 0.1109, 0.1197, 0.1305, 0.1423, 0.1541, 0.1679, 0.1835, 0.1995, 0.2181, 0.2371, 0.2578, 0.2817, 0.3061, 0.3339, 0.3661, 0.3989, 0.4348, 0.4714]

    elif(instrument == 'het'):

        if(data_type == 'll'):

            channels = [0,1,2,3]

            for i in channels:
                df_electron_fluxes = df_electron_fluxes.rename(columns={'Ele_Flux_{}'.format(i):'Electron_Flux_{}'.format(i)})
                df_electron_uncertainties = df_electron_uncertainties.rename(columns={'Ele_Flux_Sigma_{}'.format(i):'Electron_Uncertainty_{}'.format(i)})

            e_low = [0.4533, 1.0530, 2.4010, 5.9930]
            e_high = [1.0380, 2.4010, 5.9930, 18.8300]

        elif(data_type == 'l2'):

            channels = [0,1,2,3]

            e_low = [0.4533, 1.0530, 2.4010, 5.9930]
            e_high = [1.0380, 2.4010, 5.9930, 18.8300]

    elif(instrument == 'step'):

        if(data_type == 'l2'):

            channels = range(0,48)

            step_data = make_step_electron_flux(df_electrons, mask_conta=masking)
            print(step_data[0].index)
            
            df_electron_fluxes = step_data[0][plotstart:plotend]
            df_electron_uncertainties = step_data[1][plotstart:plotend]

            e_low = step_data[2]
            e_high = step_data[3]

        # Cleans up negative flux values in STEP data.
        df_electron_fluxes[df_electron_fluxes<0] = np.NaN

    if(averaging_mode == 'mean'):

        if(instrument=='ept'):

            df_proton_fluxes =df_proton_fluxes.resample('{}min'.format(averaging)).mean()
            df_proton_uncertainties = df_proton_uncertainties.resample('{}min'.format(averaging)).apply(average_flux_error)

        df_electron_fluxes = df_electron_fluxes.resample('{}min'.format(averaging)).mean()
        df_electron_uncertainties = df_electron_uncertainties.resample('{}min'.format(averaging)).apply(average_flux_error)

    # The rolling window might be broken, but it's not ever used.
    elif(averaging_mode == 'rolling_window'):

        df_electron_fluxes = df_electron_fluxes.rolling(window=averaging, min_periods=1).mean()


    if(ion_conta_corr and (instrument == 'ept')):

        ion_cont_corr_matrix = np.loadtxt('EPT_ion_contamination_flux_paco.dat')
        Electron_Flux_cont = np.zeros(np.shape(df_electron_fluxes))
        Electron_Uncertainty_cont = np.zeros(np.shape(df_electron_uncertainties))

        for tt in range(len(df_electron_fluxes)):

            # Electron_Flux_cont[tt,:] = np.sum(ion_cont_corr_matrix * df_protons.Ion_Flux.values[tt,:], axis=1)
            Electron_Flux_cont[tt, :] = np.matmul(ion_cont_corr_matrix, df_proton_fluxes.values[tt, :])
            Electron_Uncertainty_cont[tt, :] = np.sqrt(np.matmul(ion_cont_corr_matrix**2, df_proton_uncertainties.values[tt, :]**2 ))

        df_electron_fluxes = df_electron_fluxes - Electron_Flux_cont
        df_electron_uncertainties = np.sqrt(df_electron_uncertainties**2 + Electron_Uncertainty_cont**2 )
    
    if(instrument=='ept'):

        ion_string = 'Ion_contamination_correction'

    elif(instrument=='step'):

        ion_string = 'Ion_masking'

    elif(instrument=='het'):

        ion_string = ''

    # Main information dataframe containing most of the required data.
    df_info = pd.DataFrame({'Plot_period':[], 'Search_period':[], 'Bg_period':[], 'Averaging':[], '{}'.format(ion_string):[], 'Energy_channel':[], 'Primary_energy':[], 'Energy_error_low':[], 'Energy_error_high':[], 'Peak_timestamp':[], 'Flux_peak':[], 'Peak_significance':[], 'Peak_electron_uncertainty':[], 'Background_flux':[],'Bg_electron_uncertainty':[], 'Bg_subtracted_peak':[], 'Backsub_peak_uncertainty':[]})

    # Adds basic metadata to main info df.
    df_info['Plot_period'] = [plotstart]+[plotend]+['']*(len(channels)-2)
    df_info['Search_period'] = [searchstart]+[searchend]+['']*(len(channels)-2)
    df_info['Bg_period'] = [bgstart]+[bgend]+['']*(len(channels)-2)

    if(instrument=='ept'):

        df_info['Ion_contamination_correction'] = [ion_conta_corr]+['']*(len(channels)-1)

    elif(instrument=='step'):

        df_info['Ion_masking'] = [masking]+['']*(len(channels)-1)

    if(averaging_mode == 'none'):

        df_info['Averaging'] = ['No averaging']+['']*(len(channels)-1)

    elif(averaging_mode == 'rolling_window'):

        df_info['Averaging'] = ['Rolling window', 'Window size = ' + str(averaging)] + ['']*(len(channels)-2)

    elif(averaging_mode == 'mean'):

        df_info['Averaging'] = ['Mean', 'Resampled to ' + str(averaging) + 'min'] + ['']*(len(channels)-2)

    # Energy bin primary energies; geometric mean.
    primary_energies = []

    for i in range(0,len(e_low)):

        primary_energies.append(np.sqrt(e_low[i]*e_high[i]))

    primary_energies_channels = []

    for energy in channels:

        primary_energies_channels.append(primary_energies[energy])

    df_info['Primary_energy'] = primary_energies_channels

    # Next blocks of code calculate information from data and append them to main info df.
    list_bg_fluxes = []
    list_flux_peaks = []
    list_peak_timestamps = []
    list_bg_subtracted_peaks = []
    list_peak_electron_uncertainties = []
    list_average_bg_uncertainties = []
    list_bg_std = []
    list_peak_significance = []

    for channel in channels:

        bg_flux = df_electron_fluxes['Electron_Flux_{}'.format(channel)][bgstart:bgend].mean(skipna=True)
        list_bg_fluxes.append(bg_flux)

        flux_peak = df_electron_fluxes['Electron_Flux_{}'.format(channel)][searchstart:searchend].max()
        list_flux_peaks.append(flux_peak)

        peak_timestamp = df_electron_fluxes['Electron_Flux_{}'.format(channel)][searchstart:searchend].idxmax()
        list_peak_timestamps.append(peak_timestamp)

        # First finding the index location of the peak timestamp in uncertainty dataframe and the getting value of that index location.
        timestamp_loc = df_electron_uncertainties['Electron_Uncertainty_{}'.format(channel)].index.get_loc(peak_timestamp, method='nearest')
        peak_electron_uncertainty = df_electron_uncertainties['Electron_Uncertainty_{}'.format(channel)].iloc[timestamp_loc]

        #peak_electron_uncertainty = df_electron_uncertainties['Electron_Uncertainty_{}'.format(channel)][peak_timestamp]
        list_peak_electron_uncertainties.append(peak_electron_uncertainty)

        average_bg_uncertainty = np.sqrt((df_electron_uncertainties['Electron_Uncertainty_{}'.format(channel)]
                                          [bgstart:bgend]**2).sum(axis=0))/len(df_electron_uncertainties['Electron_Uncertainty_{}'.format(channel)][bgstart:bgend])
        list_average_bg_uncertainties.append(average_bg_uncertainty)

        bg_std = df_electron_fluxes['Electron_Flux_{}'.format(channel)][bgstart:bgend].std()
        list_bg_std.append(bg_std)
    
    for i in range(0,len(list_flux_peaks)):

        list_bg_subtracted_peaks.append(list_flux_peaks[i]-list_bg_fluxes[i])
        list_peak_significance.append(list_bg_subtracted_peaks[i]/list_bg_std[i])
        #sometimes the background can be higher than the peak to need to delete those values (set to nan)
        if list_bg_subtracted_peaks[i]<list_bg_fluxes[i]:
             list_peak_significance[i] = -1


    
    df_info['Energy_channel'] = channels
    df_info['Background_flux'] = list_bg_fluxes
    df_info['Flux_peak'] = list_flux_peaks
    df_info['Peak_timestamp'] = list_peak_timestamps
    df_info['Bg_subtracted_peak'] = list_bg_subtracted_peaks
    df_info['Peak_electron_uncertainty'] = list_peak_electron_uncertainties
    df_info['Bg_electron_uncertainty'] = list_average_bg_uncertainties
    df_info['Peak_significance'] = list_peak_significance
    df_info['Backsub_peak_uncertainty'] = np.sqrt(df_info['Peak_electron_uncertainty']**2 + df_info['Bg_electron_uncertainty']**2)
    


    # Calculates energy errors for spectrum plot.
    energy_error_low = []
    energy_error_high = []

    for i in range(0,len(primary_energies)):

        energy_error_low.append(primary_energies[i]-e_low[i])
        energy_error_high.append(e_high[i]-primary_energies[i])

    energy_error_low_channels = []
    energy_error_high_channels = []

    for i in channels:

        energy_error_low_channels.append(energy_error_low[i])
        energy_error_high_channels.append(energy_error_high[i])

    df_info['Energy_error_low'] = energy_error_low_channels
    df_info['Energy_error_high'] = energy_error_high_channels

    return df_electron_fluxes, df_info, [searchstart, searchend], [e_low, e_high], [instrument, data_type]

# Workaround for STEP data, there's probably a better way in Python to handle this.
def extract_step_data(df_particles, plotstart, plotend, searchstart, searchend, bgstart, bgend, instrument = 'step', data_type = 'l2', averaging_mode='none', averaging=2, masking=False, ion_conta_corr=False):

    return extract_data(df_particles, df_particles, plotstart, plotend, searchstart, searchend, bgstart, bgend, instrument = instrument, data_type = data_type, averaging_mode=averaging_mode, averaging=averaging, masking=masking, ion_conta_corr=ion_conta_corr)

def make_step_electron_flux(stepdata, mask_conta=True):
    '''
    here we use the calibration factors from Paco (Alcala) to calculate the electron flux out of the (integral - magnet) fluxes (we now use level2 data files to get these)
    we also check if the integral counts are sufficiently higher than the magnet counts so that we can really assume it's electrons (otherwise we mask the output arrays)
    As suggested by Alex Kollhoff & Berger use a 5 sigma threshold:
    C_INT >> C_MAG:
    C_INT - C_MAG > 5*sqrt(C_INT)
    Alex: die count rates und fuer die uebrigen Zeiten gebe ich ein oberes Limit des Elektronenflusses an, das sich nach 5*sqrt(C_INT) /(E_f - E_i) /G_e berechnet.
    '''
    # calculate electron flux from F_INT - F_MAG:
    colnames = ["ch_num", "E_low", "E_hi", "factors"]
    paco = pd.read_csv('step_electrons_calibration.csv', names=colnames, skiprows=1)
    paco.E_low = round(paco.E_low/1000, 5)
    paco.E_hi = round(paco.E_hi/1000, 5)

    F_INT = stepdata['Integral_Flux']
    F_MAG = stepdata['Magnet_Flux']
    step_flux =  (F_INT - F_MAG) * paco.factors.values
    U_INT = stepdata['Integral_Uncertainty']
    U_MAG = stepdata['Magnet_Uncertainty']
    # from Paco:
    # Ele_Uncertainty = k * sqrt(Integral_Uncertainty^2 + Magnet_Uncertainty^2)
    step_unc = np.sqrt(U_INT*2 + U_MAG*2) * paco.factors.values
    param_list = ['Electron_Flux', 'Electron_Uncertainty']

    if mask_conta:

        C_INT = stepdata['Integral_Rate']
        C_MAG = stepdata['Magnet_Rate']
        clean = (C_INT - C_MAG) > 5*np.sqrt(C_INT)
        step_flux = step_flux.mask(clean)
        step_unc = step_unc.mask(clean)
    step_data = pd.concat([step_flux, step_unc], axis=1, keys=param_list)

    df_electron_fluxes = step_data['Electron_Flux']
    df_electron_uncertainties = step_data['Electron_Uncertainty']

    for channel in df_electron_fluxes:

        df_electron_fluxes = df_electron_fluxes.rename(columns={channel:'Electron_Flux_{}'.format(channel)})

    for channel in df_electron_uncertainties:

        df_electron_uncertainties = df_electron_uncertainties.rename(columns={channel:'Electron_Uncertainty_{}'.format(channel)})

    return df_electron_fluxes, df_electron_uncertainties, paco.E_low, paco.E_hi

def average_flux_error(flux_err: pd.DataFrame) -> pd.Series:

    return np.sqrt((flux_err ** 2).sum(axis=0)) / len(flux_err.values)

def plot_channels(args, bg_subtraction=False, savefig=False, sigma = 3, path='', key=''):
    peak_sig = args[1]['Peak_significance']
    hours = mdates.HourLocator(interval = 1)
    df_electron_fluxes = args[0]
    df_info = args[1]
    search_area = args[2]
    energy_bin = args[3]
    instrument = args[4][0]
    data_type = args[4][1]

    title_string = instrument.upper() + ', ' + data_type.upper() + ', ' + str(df_info['Plot_period'][0][:-5])
    filename = 'channels-' + str(df_info['Plot_period'][0][:-5]) + '-' + instrument.upper() + '-' + data_type.upper() 
    
    if(df_info['Averaging'][0]=='Mean'):

        title_string = title_string + ', ' + df_info['Averaging'][1].split()[2] + ' averaging'
        filename = filename + '-' + df_info['Averaging'][1].split()[2] + '_averaging'

    elif(df_info['Averaging'][0]=='No averaging'):

        title_string = title_string + ', no averaging'
        filename = filename + '-no_averaging'

    if(bg_subtraction):
        
       title_string = title_string + ', bg subtraction on'
       filename = filename + '-bg_subtr'

    else:

        title_string = title_string + ', bg subtraction off'
    
    if(instrument == 'ept'):
        
        if(df_info['Ion_contamination_correction'][0]):

            title_string = title_string + ', ion correction on'
            filename = filename + '-ion_corr'

        elif(df_info['Ion_contamination_correction'][0]==False):

            title_string = title_string + ', ion correction off'

    # If background subtraction is enabled, subtracts bg_flux from all observations. If flux value is negative, changes it to NaN.
    if(bg_subtraction == False):
        pass
    elif(bg_subtraction == True):
        df_electron_fluxes = df_electron_fluxes.sub(df_info['Background_flux'].values, axis=1)
        df_electron_fluxes[df_electron_fluxes<0] = np.NaN

    # Plotting part.
    # Initialized the main figure.
    fig = plt.figure()
    plt.xticks([])
    plt.yticks([])
    plt.ylabel("Flux \n [1/s cm$^2$ sr MeV]", labelpad=40)
    plt.xlabel("Time", labelpad=45)
    plt.title(title_string)

    # Loop through selected energy channels and create a subplot for each.
    n=1
    for channel in df_info['Energy_channel']:

        ax = fig.add_subplot(len(df_info['Energy_channel']),1,n)
        ax = df_electron_fluxes['Electron_Flux_{}'.format(channel)].plot(logy=True, figsize=(20,25), color='red', drawstyle='steps-mid')

        plt.text(0.025,0.7, str(energy_bin[0][channel]) + " - " + str(energy_bin[1][channel]) + " MeV", transform=ax.transAxes, size=13)

        # Search area vertical lines.
        ax.axvline(search_area[0], color='black')
        ax.axvline(search_area[1], color='black')

        # Peak vertical line.
        if peak_sig[n-1]<sigma:
            ax.axvline(df_info['Peak_timestamp'][n-1], color='gray')
        if peak_sig[n-1]>sigma:
            ax.axvline(df_info['Peak_timestamp'][n-1], color='green')

        # Background measurement area.
        ax.axvspan(df_info['Bg_period'][0], df_info['Bg_period'][1], color='gray', alpha=0.25)

        ax.get_xaxis().set_visible(False)

        if(n == len(df_info['Energy_channel'])):

            ax.get_xaxis().set_visible(True)

        plt.xlabel("")
        #ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d\n%H:%M"))
        #ax.xaxis.set_minor_locator(hours)

        n+=1

    # Saves figure, if enabled.
    if(path[len(path)-1] != '/'):

        path = path + '/'

    if(savefig):

        plt.savefig(path + filename + str(key) +'.jpg', bbox_inches='tight')

    plt.show()

# This plot_check function is not finished, but it does produce cool rainbow colored plots.
def plot_check(args, bg_subtraction=False, savefig=False, key=''):

    hours = mdates.HourLocator(interval = 1)
    df_electron_fluxes = args[0]
    df_info = args[1]
    search_area = args[2]
    energy_bin = args[3]
    instrument = args[4][0]
    data_type = args[4][1]

    fig = plt.figure()
    colors = iter(plt.cm.jet(np.linspace(0, 1, len(df_info['Energy_channel']))))

    #for channel in df_info['Energy_channel']:
    #    ax = df_electron_fluxes['Electron_Flux_{}'.format(channel)].plot(logy=True, figsize=(20,25), color='red', drawstyle='steps-mid')

    for channel in df_info['Energy_channel']:

        col = next(colors)
        ax = df_electron_fluxes['Electron_Flux_{}'.format(channel)].plot(logy=True, figsize=(13,10), color=col, drawstyle='steps-mid')

    plt.show()

def plot_spectrum(args, bg_subtraction=True, savefig=False, path='', key=''):

    df_info = args[1]
    instrument = args[4][0]
    data_type = args[4][1]
    
    title_string = instrument.upper() + ', ' + data_type.upper() + ', ' + str(df_info['Plot_period'][0][:-5])
    filename = 'spectrum-' + str(df_info['Plot_period'][0][:-5]) + '-' + instrument.upper() + '-' + data_type.upper() 
    
    if(df_info['Averaging'][0]=='Mean'):

        title_string = title_string + ', ' + df_info['Averaging'][1].split()[2] + ' averaging'
        filename = filename + '-' + df_info['Averaging'][1].split()[2] + '_averaging'

    elif(df_info['Averaging'][0]=='No averaging'):

        title_string = title_string + ', no averaging'
        filename = filename + '-no_averaging'

    if(bg_subtraction):
        
       title_string = title_string + ', bg subtraction on'
       filename = filename + '-bg_subtr'

    else:

        title_string = title_string + ', bg subtraction off'
    
    if(instrument == 'ept'):

        if(df_info['Ion_contamination_correction'][0] and instrument=='ept'):

            title_string = title_string + ', ion correction on'
            filename = filename + '-ion_corr'

        elif(df_info['Ion_contamination_correction'][0]==False):

            title_string = title_string + ', ion correction off'



    # Plots either the background subtracted or raw flux peaks depending on choice.
    if(bg_subtraction):
        ax = df_info.plot.scatter(x='Primary_energy', y='Bg_subtracted_peak', c='red', label='Flux peaks', figsize=(13,10))
        ax.errorbar(x=df_info['Primary_energy'], y=df_info['Bg_subtracted_peak'], yerr=df_info['Backsub_peak_uncertainty'],
                    xerr=[df_info['Energy_error_low'], df_info['Energy_error_high']], fmt='.', ecolor='red', alpha=0.5)
    elif(bg_subtraction == False):

        ax = df_info.plot.scatter(x='Primary_energy', y='Flux_peak', c='red', label='Flux peaks', figsize=(13,10))
        ax.errorbar(x=df_info['Primary_energy'], y=df_info['Flux_peak'], yerr=df_info['Peak_electron_uncertainty'],
                    xerr=[df_info['Energy_error_low'], df_info['Energy_error_high']], fmt='.', ecolor='red', alpha=0.5)
    
    # Plots background flux and background errorbars in same scatterplot.
    df_info.plot(kind='scatter', x='Primary_energy', y='Background_flux', c='red', alpha=0.25, ax=ax, label='Background flux')
    ax.errorbar(x=df_info['Primary_energy'], y=df_info['Background_flux'], yerr=df_info['Bg_electron_uncertainty'], xerr=[df_info['Energy_error_low'],df_info['Energy_error_high']],
                fmt='.', ecolor='red', alpha=0.15)

    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlabel('Energy [MeV]', size=20)
    ax.set_ylabel('Flux \n [1/s cm$^2$ sr MeV]', size=20)
    plt.tick_params(axis='x', which='minor', labelsize=16)
    ax.xaxis.set_minor_formatter(FormatStrFormatter("%.2f"))
    #plt.tick_params(axis='y', which='minor')
    #ax.yaxis.set_minor_formatter(FormatStrFormatter("%.0f"))
    plt.legend(prop={'size': 18})
    plt.xticks(size=16)
    plt.yticks(size=16)
    plt.grid()
    plt.title(title_string)

    for label in ax.xaxis.get_ticklabels(which='minor')[1::2]:

        label.set_visible(False)
    
    if(path[len(path)-1] != '/'):

        path = path + '/'

    if(savefig):

        plt.savefig(path + filename + str(key) +'.jpg', dpi=300, bbox_inches='tight')

    plt.show()

def write_to_csv(args, path='', key=''):

    df_info = args[1]
    instrument = args[4][0]
    data_type = args[4][1]
    
    filename = 'electron_data-' + str(df_info['Plot_period'][0][:-5]) + '-' + instrument.upper() + '-' + data_type.upper()

    if(df_info['Averaging'][0] == 'Mean'):
        
        filename = filename + '-' + df_info['Averaging'][1].split()[2] + '_averaging'

    elif(df_info['Averaging'][0] == 'No averaging'):

        filename = filename + '-no_averaging'

    if(instrument == 'ept'):

        if(df_info['Ion_contamination_correction'][0]):

            filename = filename + '-ion_corr'

    df_info.to_csv(path + filename + str(key) + '.csv', index=False)

# This acc_flux function is not really finished, just something I put together quickly.
def acc_flux(args, time=[]):

    df_electron_fluxes = args[0]
    df_info = args[1]

    # If no timeframe specified, use search area.
    if(time==[]):

        time = args[2]

    # Calculates average fluxes for each enery channel from given timeframe and appends to list.
    list_flux_averages = []

    for channel in df_info['Energy_channel']:

        list_flux_averages.append(df_electron_fluxes['Electron_Flux_{}'.format(channel)][time[0]:time[1]].mean())

    df_acc = pd.DataFrame({'Primary_energy':[], 'Acc_flux':[]})
    df_acc['Primary_energy'] = df_info['Primary_energy']
    df_acc['Acc_flux'] = list_flux_averages

    ax = df_acc.plot(kind='scatter', x='Primary_energy', y='Acc_flux', logy=True, logx=True, color='green', figsize=(13,10))
