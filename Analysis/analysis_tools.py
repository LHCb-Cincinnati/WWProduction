'''
'''

def calculate_weights(cross_section, array):
    scale_factor = bool(cross_section)*(cross_section / len(array) - 1) + 1
    weights_array = [scale_factor] * len(array)
    return(weights_array)

def create_folder_path(file_name, test_mode_flag):
    import os

    if test_mode_flag:
        folder_name = 'Test'
    else:
        folder_name = file_name[:-5]
        folder_name = folder_name.split('/')[-1]
    path = os.environ['HOME'] + '/WWProduction/Data/Figures/' + folder_name
    if not os.path.exists(path):
        os.mkdir(path)
    return(path)

def calculate_hist_stats(hist, bins):
    ''' Calculate count, mean, and variance for a numpy histogram.

    Calculate count, mean, and variance of the given histogram using the
    histogram, not the data in the histogram.

    Args:
    hist (np.array): A numpy histogram that represents some dataset.
    bins (np.array): The bin edge placements for hist.

    Returns:
    hist_count (float): The number of samples in hist.
    hist_mean (float): The average value of hist.
    hist_variance (float): The variance of the samples in hist.

    '''
    import numpy as np

    hist_count = np.sum(hist)
    hist_mids = 0.5*(bins[1:] + bins[:-1])
    hist_mean = np.average(hist_mids, weights=(hist/hist_count))
    hist_var = np.sqrt(np.average((hist_mids - hist_mean)**2,
                       weights=(hist/hist_count)))
    return(hist_count, hist_mean, hist_var)

def create_hist(array, title, yscale='linear', **kwargs):
    ''' Create a 1D histogram from a numpy array and save it.

    Create a 1D histogram from a numpy array, and save it to
    a file.  The file name will be derived from the title of the histogram.
    The recommended usage of this function is to specify more keyword
    arguments than is required.
    Ex: create__hist(array1, 'Histogram 1', yscale='log', bins=50,
                    range=(0,100))

    Args:
    array (np.array): A numpy array of the data to be histogrammed.
    title (str): The title of the new histogram.
    yscale (str): The type of scale used for the yaxis of this histogram.
        Should be either 'linear' or 'log'.
    **kwargs:  Any additional keyword arguments are fed into the matplotlib 
        hist function.

    Returns:
    None

    '''
    import matplotlib.pyplot as plt

    fig, axs = plt.subplots()
    plt.subplots_adjust(top=0.85)
    hist, bins, patches = axs.hist(array, **kwargs)
    hist_count, hist_mean, hist_var = calculate_hist_stats(hist, bins)
    plt.yscale(yscale)
    plt.title(title)
    fig_string = (f"Statistics:\n"
                  f"Count: {hist_count:.2f}\n"
                  f"Mean:  {hist_mean:.2f}\n"
                  f"Sigma: {hist_var:.2f}")

    axs.text(0.8, 1.02, fig_string, transform=axs.transAxes,
              bbox=dict(facecolor='none', edgecolor='0.7', pad=3.0))

    # Slightly fancy to remove whitespace
    save_str = ''.join(title.split())
    plt.savefig(save_str + '.png')

def create_stacked_hist(array_list ,title, yscale='linear', **kwargs):
    ''' Create a 1D stacked histogram from several numpy arrays and save it.

    Create a 1D stacked histogram from several numpy arrays, and saves it to
    a file.  The file name will be derived from the title of the histogram.
    The recommended usage of this function is to specify more keyword
    arguments than is required.
    Ex: create_stacked_hist((array1, array2), 'Histogram 1', yscale='log',
                            bins=50, range=(0,100),
                            label=['array1 label', 'array2 label'])

    Args:
    array_list (list[np.array]): A list of numpy arrays full of the data to be
        histogrammed.
    title (str): The title of the new histogram.
    yscale (str): The type of scale used for the yaxis of this histogram.
        Should be either 'linear' or 'log'.
    **kwargs:  Any additional keyword arguments are fed into the matplotlib 
        hist function.

    Returns:
    None

    '''
    import matplotlib.pyplot as plt

    fig, axs = plt.subplots()
    plt.subplots_adjust(top=0.85)
    hist, bins, patches = axs.hist(array_list, stacked=True, **kwargs)
    plt.yscale(yscale)
    plt.title(title)
    #plt.legend()

    #Slightly fancy to remove whitespace
    save_str = ''.join(title.split())
    plt.savefig(save_str + '.png')