# =======================================================================================
##
# For usage: $python FatesGridComp.py -h
#
#
#
# This script is intended to diagnose the output of gridded CLM/ELM/FATES
# simulations, for rapid visualization.
#
# Options include 1) the ability to do a regression against another set of files (base).
#		  3) plotting (most of the analysis right now is visual, not much
#			       really happens right now without plotting)
#
#  UNcomment the imported "code" library and add the following call
#  in the code where you want to add a debug stop:
#	 code.interact(local=locals())
#
#
#
# =======================================================================================

import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import sys
import getopt
import code  # For development: code.interact(local=locals())
import datetime
from scipy.io import netcdf
from mpl_toolkits.basemap import Basemap
from FatesMapFunctions import map_type, map_plot_type, map_plot_title_type, PlotMaps, GetNCList
from FatesMapFunctions import pdf_type,png_type

# Some constants
g_to_Mg = 1.0e-6
m2_to_ha = 1.0e-4

ylgn_seq_cmap=mpl.cm.get_cmap('YlGn')
rdbu_div_cmap=mpl.cm.get_cmap('RdBu')

# Month string names
moname  = ['Jan','Feb','Mar','Apr','May','Jun', \
           'Jul','Aug','Sep','Oct','Nov','Dec']

def usage():
    print('')
    print('=======================================================================')
    print('')
    print(' python FatesGridComp.py -h --plotmode --regressmode')
    print('				--test-hist-file=<path> --base-hist-file=<path>')
    print('				--test-name=<text> --base-name=<text>')
    print('')
    print('	 This script is intended to diagnose the output of one or two gridded')
    print('	 runs, as a rapid pass/fail visual analysis of spatial ecosystem patterns')
    print('	 that have emerged over time.')
    print('')
    print('')
    print(' -h --help ')
    print('	    print this help message')
    print('')
    print(' --regressmode')
    print('	    [Optional] logical switch, turns on regression tests')
    print('	    against a baseline. Requires user to also set --base-rest-pref')
    print('	    default is False')
    print('')
    print(' --eval-id=<id-string>')
    print('	    a string that gives a name, or some id-tag associated with the')
    print('	    evaluation being conducted. This will be used in output file naming.')
    print('	    Any spaces detected in string will be trimmed.')
    print('')
    print(' --save-pref=<path>')
    print('     Optional.  If this path exists, instead of generating the plot')
    print('     in the window, it will save to file. The file name will have this')
    print('     prefix, with the eval-id test string appended, and then a counter')
    print('     appended to that')
    print('')
    print(' --test-hist-pref=<path>')
    print('	    the full path to the test history folder and file prefix')
    print('	    version of output')
    print(' --base-hist-pref=<path>')
    print('	    the full path to the base history folder and file prefix')
    print('	    version of output')
    print('')
    print(' --test-name=<text>')
    print('	    [Optional] a short descriptor for the test case that will be used')
    print('	    for labeling plots. The default for the test case is "test".')
    print('')
    print(' --base-name=<text>')
    print('	    [Optional] a short descriptor for the base case that will be used')
    print('	    for labeling plots. The default for the base case is "base".')
    print('')
    print('')
    print('=======================================================================')

    # ========================================================================================


def interp_args(argv):

    argv.pop(0)  # The script itself is the first argument, forget it

    ## Binary flag that turns on and off regression tests against a baseline run
    regressionmode = False
    ## history file from the test simulation
    test_h_pref = ''
    ## history file from the base simulation
    base_h_pref = ''
    ## Name of the evaluation being performed, this is non-optional
    eval_id = ''
    ## Name for plot labeling of the test case
    test_name = 'test'
    ## Name for plot labeling of the base case
    base_name = 'base'
    ## Prefix for saving plots
    save_pref = ''

    try:
        opts, args = getopt.getopt(argv, 'h',["help","regressmode", \
                                              "eval-id=","save-pref=", \
                                              "test-hist-pref=","base-hist-pref=", \
                                              "test-name=","base-name="])

    except getopt.GetoptError as err:
        print('Argument error, see usage')
        usage()
        sys.exit(2)

    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit(0)
        elif o in ("--regressmode"):
            regressionmode = True
        elif o in ("--eval-id"):
            eval_id = a
        elif o in ("--save-pref"):
            save_pref = a
        elif o in ("--test-hist-pref"):
            test_h_pref = a
        elif o in ("--base-hist-pref"):
            base_h_pref = a
        elif o in ("--test-name"):
            test_name = a
        elif o in ("--base-name"):
            base_name = a
        else:
            assert False, "unhandled option"

    if(regressionmode):
        print('Regression Testing is ON')
        if(base_h_file==''):
            print('In a regression style comparison, you must specify a')
            print('path to baseline history files. See usage:')
            usage()
            sys.exit(2)
        else:
            print('Regression Testing is OFF')

    if(test_h_pref==''):
        print('A path to history files is required input, see usage:')
        usage()
        sys.exit(2)

    if(eval_id==''):
        print('You must provide a name/id for this evaluation, see --eval-id:')
        usage()
        sys.exit(2)

    # Remove Whitespace in eval_id string
    eval_id.replace(" ","")

    return (regressionmode, eval_id, test_h_pref, base_h_pref, test_name, base_name,save_pref)


# ========================================================================================
# ========================================================================================
#					 Main
# ========================================================================================
# ========================================================================================

def main(argv):

    # -----------------------------------------------------------------------------------
    # Interpret the arguments to the script
    # -----------------------------------------------------------------------------------

    regressionmode, eval_id, test_h_pref, base_h_pref, \
        test_name, base_name, save_pref = interp_args(argv)

    # If a file prefix for output images was provided, it is assumed that
    # figures should be
    if(len(save_pref)>0):
        do_save = True
        d = datetime.datetime.today()
        datestr='c{}'.format(d.strftime("%Y%m%d-%H-%M"))
        outfile_name = "{}{}_{}.pdf".format(save_pref,eval_id,datestr)
        pdf = PdfPages(outfile_name)
        print('Saving images to: {}'.format(outfile_name))
    else:
        outfile_name = ''
        do_save = False
        pdf = []

    #    code.interact(local=locals())

    # -----------------------------------------------------------------------------------
    # Get a list of files from the test and the base, and load up
    # some coordinate data from the 1st file
    # -----------------------------------------------------------------------------------
    test_h0_list = GetNCList(test_h_pref,'h0')

    # This string is used for plot titles and text
    delta_str = 'Delta ({})-({})'.format(test_name.strip(),base_name.strip())

    # Load up the coordinate data
    fp_coords = netcdf.netcdf_file(test_h0_list[0], 'r', mmap=False)
    latvec_in = fp_coords.variables['lat'].data;
    lonvec_in = fp_coords.variables['lon'].data;

    # Load up the list of base files if using regression mode
    # Test to see if the number of files and the coordinates
    # are consistent with the test list
    if(regressionmode):
        base_h0_list = GetNCList(base_h_pref,'h0')

        if(len(test_h0_list) != len(base_h0_list)):
            print('Number of history files in the test and base are different')
            print('num test: {}'.format(len(test_h0_list)))
            print('num base: {}'.format(len(base_h0_list)))
            sys.exit(2)

        # Check to see if coordinates are the same
        fp_base_coords = netcdf.netcdf_file(test_h0_list[0], 'r', mmap=False)
        latvec_base = fp_base_coords.variables['lat'].data;
        if(abs(np.sum(latvec_in)-np.sum(latvec_base))>0.1):
            print('The latitudes of the test and base dont match')
            sys.exit(2)
        fp_base_coords.close()

    # -----------------------------------------------------------------------------------
    # Go clean up the coordinate array, do things like cut off the top
    # indices and create vertices for plotting.
    # Change coordinate system and create a re-order index array
    # -----------------------------------------------------------------------------------

    posids = np.where(lonvec_in>180.0)
    lonvec_in[posids] = -360.0+lonvec_in[posids]
    sort_ids=np.argsort(lonvec_in)
    lonvec_in=lonvec_in[sort_ids]

    # We will crop off the top and bottom of the f45, as those
    # points are at 90N and 90S? The coordinates we are creating
    # are going to be the 4 corners of each cell, so we end up having
    # one extra entry for each. Thus the longitudes end up with a size of
    # +1, and the latitude end up with a size of -2+1=-1

    nlat=latvec_in.size-1
    nlon=lonvec_in.size+1

    latvec = np.empty(nlat)
    lonvec = np.empty(nlon)

    dlon = lonvec_in[2]-lonvec_in[1]
    dlat = latvec_in[2]-latvec_in[1]
    lonvec[0] = np.maximum(-180.0,lonvec_in[0]-0.5*dlon)
    lonvec[1:-1] = 0.5*(lonvec_in[1:]+lonvec_in[0:-1])
    lonvec[-1] = np.minimum(lonvec[-2]+dlon,180.0)
    latvec = 0.5*(latvec_in[1:]+latvec_in[0:-1])

    # Create a mesh
    xv,yv = np.meshgrid(lonvec,latvec,sparse=False,indexing='xy')

    #land fraction
    landfrac = fp_coords.variables['landfrac'].data[1:-1,sort_ids]

    # Save the ocean-ids
    ocean_ids = np.where(landfrac>1.0)

    landfrac[ocean_ids]=np.nan

    # Datasets
    #
    # Required:
    # PFTbiomass
    # PFTleafbiomass


    # Use the Coordinate file to allocate space to variables of interest.

    biomass_test = np.zeros(landfrac.shape)
    tlai_test     = np.zeros(landfrac.shape)


    atime=0
    for ifile in range(0,len(test_h0_list)):

        fpdata = netcdf.netcdf_file(test_h0_list[ifile], 'r', mmap=False)

        # Load date information
        mcdate_raw = fpdata.variables['mcdate'].data

        for itime in range(0,mcdate_raw.shape[0]):

            atime=atime+1

            # Save the current date integer to a float
            # YYYYMMDD
            fdate=float(mcdate_raw[itime])

            yr   = np.floor(fdate/10000.0)
            mo   = np.floor((fdate-yr*10000.0)/100.0)
            dom  = fdate - (yr*10000.0 + mo*100.0)

            print('{} {}-{}'.format(int(yr),moname[int(mo)-1],int(dom)))

            biomass_test = biomass_test + g_to_Mg / m2_to_ha * \
                           np.transpose(fpdata.variables['ED_biomass'].data[itime,1:-1,sort_ids])

            tlai_test = tlai_test +  np.transpose(fpdata.variables['TLAI'].data[itime,1:-1,sort_ids])

            fpdata.close()

    biomass_test = biomass_test / float(atime)
    biomass_test[ocean_ids] = np.nan    # Set ocean-cells to nan

    tlai_test = tlai_test / float(atime)
    tlai_test[ocean_ids] = np.nan



    # Set up the plot objects

    proj_type    = 'robin'
    #    proj_type    = 'cyl'
    plot_title   = ''

    title_obj = map_plot_title_type(plot_title,0.25,0.55,14)
    plot_obj  = map_plot_type(xv,yv,proj_type,outfile_name,do_save,pdf_type,pdf)

    map_list = []
    map1=map_type(biomass_test, 'Total Biomass [MgC/ha]', ylgn_seq_cmap, [])
    map_list.append(map1)
    PlotMaps(plot_obj,map_list,title_obj)

    map_list = []
    map2=map_type(tlai_test, 'LAI', ylgn_seq_cmap, [0.,6.])
    map_list.append(map2)
    PlotMaps(plot_obj,map_list,title_obj)

    # This is same as pdf.close()
    if(do_save):
        plot_obj.save_obj.close()

# =======================================================================================
# This is the actual call to main

if __name__ == "__main__":
                                    main(sys.argv)
