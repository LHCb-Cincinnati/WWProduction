'''  Provides a series of functions and classes to use my physics analysis.

Provides a series of functions and classes to use my physics analysis.

To-Dos:
- Make a function that grabs all the normal kinematic variables from a ROOT
    Tree.
'''

# Imports
# Standard Library Packages
import argparse 
import pdb
import os
from collections import namedtuple
# 3rd Party Packages
import numpy as np
import matplotlib.pyplot as plt

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
        self.add_argument('-c', '--cross_section', type=float, default=False,
                            nargs='?', help='''The cross section used to
                            normalize histograms in fb.''')
        self.add_argument('-l', '--luminosity', type=float, default=False, 
                            nargs='?', help='''Luminosity used to scale the
                            cross-section.  The default units are in fb^-1.
                            This option can be defaulted to False, but will
                            draw an error if the cross-section is specified
                            without a luminosity.''')
        self.add_argument('-t', '--testing', default=True,
                            action='store_true', help='''Flag to indicate if
                            a new output folder should be created for this
                            analysis.  -t means that all plots will go in
                            folder labelled Test.''')
        self.add_argument('-p', '--production', dest='testing',
                            action='store_false', help='''Flag to indicate if
                            a new output folder should be created for this
                            analysis.  -p means that all plots will go in
                            a designated folder.''')
        self.add_argument('-o', '--output', default= "Test",
                            help=''' The name of a folder to
                            store the output in.''')
        self.add_argument('-d', '--debug',
                            action='store_true', default=False, 
                            help='''Flag to indicate if
                            the program should be run in debug mode.''')
        self.args = self.parse_args(sys_args)
        self.min_num_files = min_num_files
        self.check_args()
        print(self)

    def __str__(self):
        file_str = "Input Files: \n"
        for file in self.args.input_files:
            file_str += f"    ./{file.name} \n"
        xsection_str = f"Cross Section: {self.args.cross_section} \n"
        lum_string = f"Luminosity : {self.args.luminosity}fb^-1 \n"
        test_str = f"Testing Flag: {self.args.testing} \n"
        output_str = f"Output Folder: {self.args.output} \n"
        debug_str = f"Debug Flag: {self.args.debug}"
        return(file_str+xsection_str+lum_string+test_str+output_str+debug_str)



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
        
        self.check_min_num_files()
        self.check_xsection_luminosity_given()
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

    def check_xsection_luminosity_given(self):
        ''' Checks that the cross-section and luminosity are both given or not.

        Checks that the cross-section and luminosity are either both given or
        both not.  If not, an error is raised.

        Args:
        None

        Returns:
        None

        '''

        if (bool(self.args.cross_section) != bool(self.args.luminosity)):
            raise RuntimeError(f'''Either the luminosity and cross-section
                must be both specified or both not specified.  Their current
                values are {self.args.cross_section} and
                {self.args.luminosity}.''')
        return()

def create_folder_path(folder_name, test_mode_flag):
    ''' Creates a folder in ../dileptonww/Figures/ with the given folder_name.

    Creates a folder in ../dileptonww/Figures/ with the given folder_name.  If 
    test_mode_flag is specified the file_name defaults to 'Test'.

    Args:
    file_name (str): Name of the folder to be created.
    test_mode_flag (Bool): Boolean to describe if a test is being called.

    Returns:
        None
    '''

    if test_mode_flag:
        folder_name = 'Test'
    path = find_WW_path() + '/Figures/' + folder_name
    if not os.path.exists(path):
        os.mkdir(path)
    return(path)

def find_WW_path():
    ''' Finds the path to WWProduction the folder within the project.

    Finds the path to WWProduction the folder within the project.  Many paths
    are based upon where this folder is located.

    Args:
    None

    Returns:
    path (str): Path to WWProduciton folder.
    '''

    cwd_list = os.getcwd().split('/')
    WW_index = cwd_list.index('WWProduction')
    path = '/'.join(cwd_list[:WW_index+1])
    return(path)

def calculate_hist_stats(hist):
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

    hist_count = hist.sum().value
    hist_mids = 0.5*(hist.axes[0].edges[1:] + hist.axes[0].edges[:-1])
    hist_mean = np.average(hist_mids, weights=(hist.view().value/hist_count))
    hist_var = np.sqrt(np.average((hist_mids - hist_mean)**2,
                       weights=(hist.view().value/hist_count)))
    return(hist_count, hist_mean, hist_var)

def create_stair(hist, title, yscale='linear', luminosity=False, normalize=False, **kwargs):
    ''' Plots a 1D BH histogram and saves it.

    Plots a 1D BH histogram, and saves it to
    a file.  The file name will be derived from the title of the histogram.
    Ex: create_stair(array1, 'Histogram 1', yscale='log')

    Args:
    hist (bh.Histogram): A histogram from the boosted histogram package.
    title (str): The title of the new histogram.
    yscale (str): The type of scale used for the yaxis of this histogram.
        Should be either 'linear' or 'log'.
    luminosity (float): Adds a luminosity texbox with input luminosity in the
        text.  Default unit is fb^-1.
    **kwargs:  Any additional keyword arguments are fed into the matplotlib 
        hist function.

    Returns:
    None

    '''

    fig, axs = plt.subplots()
    plt.subplots_adjust(top=0.85)
    if normalize:
        hist_sum = hist.view().value.sum()
        plt.stairs(hist.view().value / hist_sum, edges=hist.axes[0].edges, **kwargs)
    else: 
        plt.stairs(hist.view().value, edges=hist.axes[0].edges, **kwargs)
    hist_count, hist_mean, hist_var = calculate_hist_stats(hist)
    plt.yscale(yscale)
    plt.title(title)
    fig_string = (f"Statistics:\n"
                  f"Count: {hist_count:.2f}\n"
                  f"Mean:  {hist_mean:.2f}\n"
                  f"Sigma: {hist_var:.2f}")
    axs.text(0.8, 1.02, fig_string, transform=axs.transAxes,
              bbox=dict(facecolor='none', edgecolor='0.7', pad=3.0))
    if luminosity:
        axs.text(0, 1.01, "$\mathcal{L} = " + f"{luminosity}" + "fb^{-1}$",
                transform=axs.transAxes, fontsize=12)
    # Slightly fancy to remove whitespace
    save_str = ''.join(title.split())
    plt.savefig(save_str + '.png')
    plt.close()

def create_stacked_stair(hist_list, title, label_list, yscale='linear', luminosity = False, normalize=False, **kwargs):
    '''Stacks several 1D BH histograms into a single plot.

    Stacks several 1D BH histograms into a single plot, and saves it to
    a file.  The file name will be derived from the title of the histogram.
    Ex: create_stacked_stair([hist1, hist2], 'Histogram 1', ["hist1", "hist2"],
                            yscale='log', luminosity=5.4, normalize=False)

    Args:
    hist_list list[bh.Histogram]: A list of histograms from the boosted
        histogram package.
    title (str): The title of the new histogram.
    label_list list[str]: A list of strings for the legend of the plot.  The
        number of histograms in hist_list should match the labels here.
    yscale (str): The type of scale used for the yaxis of this histogram.
        Should be either 'linear' or 'log'.
    normalize (bool): Boolean to decide if the histograms should be normalized.
        If yes, they are normalized to their own sum.
    luminosity (float): Adds a luminosity texbox with input luminosity in the
        text.  Default unit is fb^-1.
    **kwargs:  Any additional keyword arguments are fed into the matplotlib 
        stair function.

    Returns:
    None

    '''

    fig, axs = plt.subplots()
    plt.subplots_adjust(top=0.85)
    for index, hist in enumerate(hist_list):
        if normalize:
            hist_sum = hist.view().value.sum()
            plt.stairs(hist.view().value / hist_sum, edges=hist.axes[0].edges,
                        label=label_list[index], **kwargs)
        else: 
            plt.stairs(hist.view().value, edges=hist.axes[0].edges,
                        label=label_list[index], **kwargs)
    if luminosity:
        axs.text(0.8, 1.01, "$\mathcal{L} = " + f"{luminosity}" + "fb^{-1}$",
                transform=axs.transAxes, fontsize=12)
    plt.yscale(yscale)
    plt.title(title)
    plt.legend()
    # Slightly fancy to remove whitespace
    save_str = ''.join(title.split())
    plt.savefig(save_str + '.png')
    plt.close()
