'''  Provides a series of functions and classes to use in my physics analysis.

Provides a series of functions and classes to use in my physics analysis.
'''

# Imports
# Standard Library Packages
import argparse 
import pdb
import os
# 3rd Party Packages
import numpy as np
import matplotlib.pyplot as plt
# HEP Packages

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
    path = find_WW_path() + '/GenLevelStudies/Figures/' + folder_name
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

    hist_count = np.sum(hist)
    hist_mids = 0.5*(bins[1:] + bins[:-1])
    hist_mean = np.average(hist_mids, weights=(hist/hist_count))
    hist_var = np.sqrt(np.average((hist_mids - hist_mean)**2,
                       weights=(hist/hist_count)))
    return(hist_count, hist_mean, hist_var)

def create_stair(array, title, yscale='linear', luminosity=False, normalize=False, **kwargs):
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
    normalize (bool): Boolean to decide if the histograms should be normalized.
        If yes, they are normalized to their own sum.
    luminosity (float): Adds a luminosity texbox with input luminosity in the
        text.  Default unit is fb^-1.
    **kwargs:  Any additional keyword arguments are fed into the matplotlib 
        hist function.

    Returns:
    None

    '''

    fig, axs = plt.subplots()
    plt.subplots_adjust(top=0.85)
    hist, bins, patches = axs.hist(array, **kwargs)
    plt.cla()
    if normalize:
        hist_sum = np.sum(hist)
        axs.stairs(hist / hist_sum, edges=bins)
    else: 
        axs.stairs(hist, edges=bins)
    hist_count, hist_mean, hist_var = calculate_hist_stats(hist, bins)
    plt.yscale(yscale)
    plt.title(title)
    fig_string = (f"Statistics:\n"
                  f"Count: {hist_count:.2f}\n"
                  f"Mean:  {hist_mean:.2f}\n"
                  f"Sigma: {hist_var:.2f}"
                )
    axs.text(0.8, 1.02, fig_string, transform=axs.transAxes,
              bbox=dict(facecolor='none', edgecolor='0.7', pad=3.0))
    if luminosity:
        axs.text(0, 1.01, "$\mathcal{L} = " + f"{luminosity}" + "fb^{-1}$",
                transform=axs.transAxes, fontsize=12)
    # Slightly fancy to remove whitespace
    save_str = ''.join(title.split())
    plt.savefig(save_str + '.png')
    plt.close()

def create_stacked_stair(array_list, title, label_list, yscale='linear', normalize=False, **kwargs):
    ''' Create a 1D histogram from a numpy array and save it.

    Create a 1D histogram from a numpy array, and save it to
    a file.  The file name will be derived from the title of the histogram.
    The recommended usage of this function is to specify more keyword
    arguments than is required.
    Ex: create__hist(array1, 'Histogram 1', yscale='log', bins=50,
                    range=(0,100))

    Args:
    array_list (list[np.array]): A list of  numpy arrays containing the data
        to be histogrammed.
    title (str): The title of the new histogram.
    yscale (str): The type of scale used for the yaxis of this histogram.
        Should be either 'linear' or 'log'.
    **kwargs:  Any additional keyword arguments are fed into the matplotlib 
        hist function.

    Returns:
    None

    '''

    fig, axs = plt.subplots()
    plt.subplots_adjust(top=0.85)
    for index, array in enumerate(array_list):
        hist, bins = np.histogram(array, **kwargs)
        if normalize:
            hist_sum = np.sum(hist)
            plt.stairs(
                hist / hist_sum,
                edges=bins,
                label=label_list[index]
            )
        else: 
            plt.stairs(
                hist / hist_sum,
                edges=bins,
                label=label_list[index]
            )
    plt.yscale(yscale)
    plt.title(title)
    plt.legend()

    # Slightly fancy to remove whitespace
    save_str = ''.join(title.split())
    plt.savefig(save_str + '.png')
    plt.close()

def fill_array(array, event, index):
    ''' Fills a numpy array with particle information from a pythia event.

    Fills a given numpy array with information about a specific particle
    from a pythia event..

    Args:
    array (np.array): Numpy array to be filled.  Should have 12 elements.
    event (pythia.event): Pythia event log.
    index (int): Index of the specific particle to be looked at within the
        pythia event.

    Returns:
    array (np.array): Filled Numpy array with information from the specific
                        particle.
    '''

    particle = event[index]
    px = particle.px()
    py = particle.py()
    pz = particle.pz()
    p = np.sqrt(px*px + py*py + pz*pz)
    array[:] = (px, py, pz,
                particle.pT(), p, particle.eta(), particle.e(),
                particle.phi(), particle.m0(), particle.id(),
                particle.charge(), particle.status())
    return(array)

def get_target_children(target_index, event, child_index_list=[]):
    ''' Gets index of non-radiative target children in a pythia event.

    Finds the indicies of the children of a given target particle within a 
    given pythia event.  Children that are radiative decays 
    (target -> target + photon) are ignored.

    Args:
    target_index (int): Index of target particle in a pythia event.
    event (pythia.event): Pythia event log.
    child_index_list (list[int]): List of indices of the daughters of the 
        target particle.  This list is given as a optional argument for use 
        in recursive function calls and should generally not be needed by 
        the user.

    Returns:
    child_index_list (list[int]): As described above.
    '''

    target = event[target_index]
    for idaughter in target.daughterList():
        child_index_list.append(idaughter)
        if abs(event[idaughter].id()) == abs(target.id()):
            child_index_list = get_target_children(idaughter, event,
                                                child_index_list)
    return(child_index_list)

def check_mother_pid(event, particle, pid):
    ''' Checks if a particles mother has a given pid.

    Checks if a particle's mother has a given pid.  Radiative decays 
    (target -> target + photon) are ignored.

    Args:
    event (pythia.event): Pythia event log.
    particle (pythia.particle): Pythia particle to be investigated.
    pid (int): PID the mother of particles needs to be chcked against.

    Returns:
    bool: True if particle's mother has the given pid.  False if not.
    '''

    imother = particle.mother1()
    mother = event[imother]
    if mother.id() == pid:
        return(True)
    elif mother.id() == particle.id():
        return(check_mother_pid(event, mother, pid))
    else:
        return(False)