'''  Provides a series of functions and classes to use my physics analysis.

Provides a series of functions and classes to use my physics analysis.

To-Dos:
- Make a funciton that grabs all the normal kinematic variables from a ROOT
    Tree.
'''

import argparse 


def calculate_weights(array, cross_section, nsim=False):
    ''' Creates the weights array for a numpy histogram, given information
    on the particle creation process.

    There are three cases for this function.
    - First, no cross section is specified.  In this case the weights array
    should be filled with ones.
    - Two, a cross section is specified, but the number of simulated events
    is not provided.  In this case each weight should be the cross section
    divided by the number of events in the array.
    - Three, both the cross section and the number of simulated events are
    provided.  In this case, the weights should be the cross section divided
    by the number of simulated events.

    Args:
    array (np.array): Numpy array of data corresponding to the process of the
        cross-section.
    cross_section (float): The cross-section (in fb )of the process through
        which the data was created. Optional
    nsim (int): The number of events simulated on the given process type.
        Optional

    Returns:
    weights_array (np.array): An array of weights for each element of the
        input array.

    '''
    if not nsim:
        nsim = len(array)
    scale_factor = bool(cross_section)*((cross_section / nsim) - 1) + 1
    weights_array = [scale_factor] * len(array)
    return(weights_array)

def create_folder_path(file_name, test_mode_flag):
    import os

    if test_mode_flag:
        folder_name = 'Test'
    else:
        folder_name = file_name[:-5]
        folder_name = folder_name.split('/')[-1]
    path = find_WW_path() + '/Data/Figures/' + folder_name
    if not os.path.exists(path):
        os.mkdir(path)
    return(path)

class Parser(argparse.ArgumentParser):
    """ A user input parser for my analysis.

    A user input parser for my analysis.  Parses input files, cross sections,
    and testing flags.

    Attributes:
        args (dict): Dictionary of user inputs.  Derived from sys.argv
        min_num_files (int): The minimum number of input files allowed by the
            parser.

    """

    def __init__(self, sys_args, min_num_files = 1):
        ''' Initialization method of class.

        Initialization method of class.  Performs initial assignment using
        the argparse.ArgumentParser class __init__ method, parses the user
        input arguments, and checks the user input for any errors.

        Args:
        sys_args (list): Input from sys.argv.
        min_num_files (int): The minimum number of input files allowed by the
            parser.

        Returns:
        None

        '''

        argparse.ArgumentParser.__init__(self, description='''Process files
                                        and settings for analysis.''',
                                        prog='program.py',
                                        usage='%(prog)s file1.root'
                                         +' file2.root -c 35.6 37.8 -p')
        self.add_argument('input_files',type=open, nargs='+',
                        help='''The file or files to be anayzed. 
                        Input files should be root files.''')
        self.add_argument('-c', '--cross_section', type=float, default=False, nargs='*', 
                            help='''The cross section used to normalize
                            histograms in fb.  This option should either not
                            be specified or have exactly as many arguments as
                            input files.''')
        self.add_argument('-t', '--testing', default=True,
                            action='store_true', help='''Flag to indicate if
                            a new output folder should be created for this
                            analysis.  -t means that all plots will go in
                            folder labelled Test.''')
        self.add_argument('-p', '--production', help=''' Flag to indicate if
                            jet data is available in the user input file.''')
        self.add_argument('-j', '--jets',
                            action='store_true', default=False, 
                            help='''Flag to indicate if
                            a new output folder should be created for this
                            analysis.  -p means that all plots will be put
                            into a new output folder withthe same name as
                            the input file.''')
        self.args = self.parse_args(sys_args)
        self.min_num_files = min_num_files
        self.check_args()
        print(self.args)

    def check_args(self):
        ''' Checks self.args input for any potential errors after parsing.

        Checks self.args input for any potential errors after parsing.  This
        may include an incorrect number of files or cross sections given.
        This function is used to call other functions which perform
        individual checks.

        Args:
        None

        Returns:
        None

        '''
    
        self.check_cross_sections()
        self.check_min_num_files()
        return()
    
    def check_cross_sections(self):
        ''' Checks that the number of cross sections is correct.

        Checks that the number of cross sections given is equal to the
        number of files given.  If not, an error is raised.

        Args:
        None

        Returns:
        None

        '''

        if (not self.args.cross_section):
            self.args.cross_section = [False] * len(self.args.input_files)
        elif (len(self.args.cross_section) != len(self.args.input_files)):
            raise RuntimeError('''The number of given cross sections is
                                different from the number of given input 
                                files.  Please either don't specify a cross
                                section arugment or specify as many cross 
                                sections as there are input files.''')
        return()

    def check_min_num_files(self):
        ''' Checks that the number of input files is acceptable.

        Checks that the number of input files given is equal to or greater
        than the min_num_files parameter.  If not, an error is raised.

        Args:
        None

        Returns:
        None

        '''

        if (len(self.args.input_files) < self.min_num_files):
            raise RuntimeError(f'''The number of files given is less than the
            the required number of files.  Please give at least 
            {self.min_num_files} files as input.''')
        return()

def find_WW_path():
    ''' Finds the path to WWProduction the folder within the project.

    Finds the path to WWProduction the folder within the project.  Many paths
    are based upon where this folder is located.

    Args:
    None

    Returns:
    path (str): Path to WWProduciton folder.
    '''
    import os

    cwd_list = os.getcwd().split('/')
    WW_index = cwd_list.index('WWProduction')
    path = '/'.join(cwd_list[:WW_index+1])
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
    plt.legend()

    #Slightly fancy to remove whitespace
    save_str = ''.join(title.split())
    plt.savefig(save_str + '.png')

def create_2D_hist(array1, array2, **kwargs):
    import matplotlib.pyplot as plt
    import numpy as np
    fig, axs = plt.subplots()
    hist, x_edges, y_edges = np.histogram2d(array1, array2, **kwargs)
    X, Y = np.meshgrid(x_edges, y_edges)
    color_mesh = axs.pcolormesh(X, Y, hist)
    fig.colorbar(color_mesh, ax=axs)
    plt.title('Test')
    plt.xlabel('xlabel')
    fig.savefig('Test' + '.png')