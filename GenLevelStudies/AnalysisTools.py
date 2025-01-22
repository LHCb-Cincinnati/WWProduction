'''  Provides a series of functions and classes to use in my physics analysis.

Provides a series of functions and classes to use in my physics analysis.
'''

# Imports
# Standard Library Packages
import argparse 
import pdb
import os
import warnings
import logging
# 3rd Party Packages
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import boost_histogram as bh
# HEP Packages
import mplhep as hep

# General Plotting Settings
matplotlib.use('AGG') # Change matplotlib backend
# Load style sheet
plt.style.use(hep.style.LHCb2)  # or ATLAS/LHCb2

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

    Raise:
    Warning: If a histogram has no entries in it.
    '''

    hist_count = hist.sum().value
    hist_mids = 0.5*(hist.axes[0].edges[1:] + hist.axes[0].edges[:-1])
    if hist_count == 0:
        warnings.warn(f"Histogram has no entries in it:", stacklevel=3)
        return(0, 0, 0)
    else:
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

    # Load style sheet
    plt.style.use(hep.style.LHCb2)  # or ATLAS/LHCb2

    fig, axs = plt.subplots()
    plt.subplots_adjust(top=0.85)
    if normalize:
        hist_vals_array = hist.view().value / hist.view().value.sum()
        hist_errs_array = np.sqrt(hist.view().variance) / hist.view().value.sum()
    else: 
        hist_vals_array = hist.view().value
        hist_errs_array = np.sqrt(hist.view().variance)
    plt.bar(
        hist.axes[0].centers,
        hist_vals_array,
        width = hist.axes[0].widths,
        yerr = hist_errs_array,
        fill = False,
        **kwargs
        )
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
    
    logging.debug(
        "Arguments given to create_stacked_hist: \n"
        + f"Histograms: {hist_list}\n"
        + f"Title: {title}\n"
        + f"Label List: {label_list}\n"
        + f"yscale: {yscale}\n"
        + f"Luminosity: {luminosity}\n"
        + f"Normalize: {normalize}\n"
        + f"kwargs: {kwargs}")
    fig, axs = plt.subplots()
    plt.subplots_adjust(top=0.85)
    for index, hist in enumerate(hist_list):
        if normalize:
            hist_sum = hist.view().value.sum()
            axs.stairs(hist.view().value / hist_sum, edges=hist.axes[0].edges,
                        label=label_list[index], **kwargs)
        else: 
            axs.stairs(hist.view().value, edges=hist.axes[0].edges,
                        label=label_list[index], **kwargs)
    if luminosity:
        axs.text(0.8, 1.01, "$\mathcal{L} = " + f"{luminosity}" + "fb^{-1}$",
                transform=axs.transAxes, fontsize=12)
    axs.set_yscale(yscale)
    axs.set_title(title)
    hist_handles, hist_labels = axs.get_legend_handles_labels()
    axs.legend(hist_handles, hist_labels, loc = "lower center")
    # Slightly fancy to remove whitespace
    save_str = ''.join(title.split())
    fig.savefig(save_str + '.png')
    plt.close()

def plot_2dhist(hist, hist_name,
    yscale='linear', luminosity=False, normalize=False,
    **kwargs
    ):

    if normalize:
        norm_factor = np.sum(hist.view().value)
    else:
        norm_factor = 1
    colormesh = plt.pcolormesh(
        hist.axes[0].edges,
        hist.axes[1].edges,
        hist.view().value.T / norm_factor,
        **kwargs
    )
    colorbar = plt.colorbar(colormesh)
    plt.title(hist_name)
    plt.yscale(yscale)
    if luminosity:
        plt.text(0, 1.01, "$\mathcal{L} = " + f"{luminosity}" + "fb^{-1}$",
            transform=plt.transAxes, fontsize=24)
    save_str = ''.join((hist_name + "_" + yscale).split())
    plt.savefig(save_str + '.png')
    plt.close()


def get_index(flattened_index, size_tuple):
    index_list = [0] * len(size_tuple)
    for dim in range(len(size_tuple)-1):
        #  This won't work in three dimensions.
        #  I will need to divide the flattened index by the skipped
        #  size dimensions.
        index_list[dim] = int(flattened_index / np.prod(size_tuple[dim+1:]))
    index_list[-1] = int(flattened_index % size_tuple[-1])
    return(tuple(index_list))

def divide_bh_histograms(num_hist, denom_hist, binomial_error=False):
    """ Divides two bh histograms with binomial errors.

    Divides two bh histograms, assuming the same binning scheme, binwise.
    Independent or binomial errors can be used.

    Args:
        num_hist (bh.histogram): Numerator histogram.
        denom_hist (bh.histogram): Denominator histogram.
        binomial_errors (bool): Whether to use binomial errors.
    
    Returns:
        div_hist (bh.histogram): num_hist / denom_hist
    """
    div_hist = num_hist.copy()
    original_shape = div_hist.view().shape
    # div_hist.reset()
    # div_hist.view().reshape(len(div_hist.view().flatten()))
    grouped_array = zip(
        range(len(div_hist.view().flatten().value)),
        num_hist.view().flatten().value,
        denom_hist.view().flatten().value,
        num_hist.view().flatten().variance,
        denom_hist.view().flatten().variance
    )
    for (index, num_value, denom_value, num_var, denom_var) in grouped_array:
        if denom_value == 0:
            continue
            warnings.warn(f"Denominator value of zero in hist: {denom_hist}")
        div_value = num_value / denom_value
        # print(div_val_array[index])
        if binomial_error:
            div_var = (
                (1/denom_value)
                * np.sqrt(
                    num_value
                    * (1- (num_value / denom_value))
                )
            )
        else:
            div_var = (
                (num_value / denom_value)
                * np.sqrt(
                    (np.sqrt(num_var) / num_value)**2
                    + (np.sqrt(denom_var) / denom_value)**2
                )
            )
        index_list = get_index(index, original_shape)
        div_hist[index_list] = [div_value, div_var]
    # div_hist.view().reshape(original_shape)
    return(div_hist)


    # div_hist = bh.Histogram(
    #     bh.axis.Variable(num_hist.axes[0].edges),
    #     storage=bh.storage.Weight()
    # )
    # for ihist_bin in range(len(num_hist.view())):
    #     num_value = num_hist[ihist_bin].value
    #     num_var = num_hist[ihist_bin].variance
    #     denom_value = denom_hist[ihist_bin].value
    #     denom_var = denom_hist[ihist_bin].variance
    #     if denom_value == 0:
    #         continue
    #         warnings.warn(f"Denominator value of zero in hist: {denom_hist}")
    #     bin_value = num_value / denom_value
    #     if binomial_error:
    #         bin_error = (
    #             (1/denom_value)
    #             * np.sqrt(
    #                 num_value
    #                 * (1- (num_value / denom_value))
    #             )
    #         )
    #         other_form = (
    #             np.abs((((1 - 2 * num_value / denom_value) * num_var) + (num_value**2 * denom_var / denom_value**2)) / denom_value**2)
    #         )
    #         print(f"Bin: {ihist_bin}")
    #         print(f"Numerator Content: {num_hist[ihist_bin].value}")
    #         print(f"Denominator Content: {denom_hist[ihist_bin].value}")
    #         print(f"Numerator Variance: {num_hist[ihist_bin].variance}")
    #         print(f"Denominator Variance: {denom_hist[ihist_bin].variance}")
    #         print(f"Normal: {bin_error}")
    #         print(f"New: {other_form}")
    #     else:
    #         bin_error = (
    #             (num_value / denom_value)
    #             * np.sqrt(
    #                 (np.sqrt(num_var) / num_value)**2
    #                 + (np.sqrt(denom_var) / denom_value)**2
    #             )
    #         )
    #     div_hist[ihist_bin] = [bin_value, bin_error]
    return(div_hist)

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

def fill_array_big(array, event, index):
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
                particle.charge(), particle.status(),
                particle.hasVertex(),
                particle.xProd(), particle.yProd(), particle.zProd(),
                particle.tProd())
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

def check_radiative_decay(event, particle):
    ''' Checks if a particle decays radiatively.
    '''

    idaughter1 = particle.daughter1()
    idaughter2 = particle.daughter2()
    if idaughter1 == 0:
        return(False)
    elif idaughter2 == 0:
        if event[idaughter1] == particle.id():
            return(True)
        else:
            return(False)
    else:
        for index in range(idaughter2 - idaughter1):
            if event[idaughter1 + index] == particle.id():
                return(True)
        return(False)

