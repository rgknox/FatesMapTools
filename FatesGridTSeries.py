# =======================================================================================
##
# For usage: $python acre_gridcomp.py -h
#
#
#
# This script is intended to diagnose the output of gridded CLM/ELM/FATES
# simulations, for rapid visualization.
#
# Options include 1) the ability to do a regression against another set of files (base).
#                 3) plotting (most of the analysis right now is visual, not much
#                              really happens right now without plotting)
#
#  UNcomment the imported "code" library and add the following call
#  in the code where you want to add a debug stop:
#        code.interact(local=locals())
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
import code  # For development: code.interact(local=dict(globals(), **locals())) 
import time
from scipy.io import netcdf
from mpl_toolkits.basemap import Basemap
import matplotlib.animation as animation

# Some constants
g_to_Mg = 1.0e-6
m2_to_ha = 1.0e-4

ylgn_seq_cmap=mpl.cm.get_cmap('YlGn')
rdbu_div_cmap=mpl.cm.get_cmap('RdBu')

def usage():
	print('')
	print('=======================================================================')
	print('')
	print(' python acre_gridcomp.py -h --test-hist-pref=<text>')
	print('                         --eval-id=<text> --var-name=<text>')
	print('')
	print('  This script is intended to diagnose the output of one or two gridded')
	print('  runs, as a rapid pass/fail visual analysis of spatial ecosystem patterns')
	print('  that have emerged over time.')
	print('')
	print('')
	print(' -h --help ')
	print('     print this help message')
	print('')
	print(' --eval-id=<id-string>')
	print('     a string that gives a name, or some id-tag associated with the')
	print('     evaluation being conducted. This will be used in output file naming.')
	print('     Any spaces detected in string will be trimmed.')
	print('')
	print(' --test-hist-pref=<path>')
	print('     the full path to the history file folder')
	print('     version of output')
	print('')
	print(' --var-name=<NETCDF_VARIABLE>')
	print('     the name of the netcdf variable of interest')
	print('     MUST BE A 1D SITE VARIABLE (for now)')
	print('')
	print('')
	print('=======================================================================')


# ========================================================================================

## interp_args processes the arguments passed to the python script at executions.
#  The options are parsed and logical checks and comparisons are made.
# @param argv
# @return plotmode
# @return regressionmode
# @return restartmode
# @return test_r_prefix
# @return base_r_prefix
# @return test_h_prefix
# @return base_h_prefix
# @return test_name
# @return base_name

def interp_args(argv):

	argv.pop(0)  # The script itself is the first argument, forget it

	## history file from the test simulation
	test_h_pref = ''
	## Name of the evaluation being performed, this is non-optional
	eval_id = ''
	## Name of the variable of interest
	var_name = ''

	try:
		opts, args = getopt.getopt(argv, 'h',["help", \
                                              "eval-id=", \
                                              "test-hist-pref=", \
                                              "var-name="])

	except getopt.GetoptError as err:
		print('Argument error, see usage')
		usage()
		sys.exit(2)

	for o, a in opts:
		if o in ("-h", "--help"):
			print('Called help')
			usage()
			sys.exit(0)
		elif o in ("--eval-id"):
			eval_id = a
		elif o in ("--test-hist-pref"):
			test_h_pref = a
		elif o in ("--var-name"):
			var_name = a
		else:
			assert False, "unhandled option"

	if(test_h_pref==''):
		print('A path to history files is required input, see usage:')
		usage()
		sys.exit(2)

	if(eval_id==''):
		print('You must provide a name/id for this evaluation, see --eval-id:')
		usage()
		sys.exit(2)

	if(var_name==''):
		print('You must provide a variable name, see usage:')
		usage()
		sys.exit(2)

	# Remove Whitespace in eval_id string
	eval_id.replace(" ","")
    
	return (eval_id, test_h_pref, var_name )


# =======================================================================================

## get a list of files given a directory prefix and specification on the type
# @param file_prefix a string with the full or relative path to the data
# @return filelist a list of strings, each of which a file

def getnclist(file_prefix,filetype):

	from os.path import isfile, join 
	from os import listdir

	def FindChar(s, ch):
		return [i for i, char in enumerate(s) if char == ch]

	slashpos = FindChar(file_prefix,'/')

	if(len(slashpos)==0):
		pathpref=file_prefix+'/'
		filepref=''
	else:
		pathpref=file_prefix[:slashpos[-1]+1]
		filepref=file_prefix[slashpos[-1]+1:]

	filelist = sorted([pathpref+f for f in listdir(pathpref) \
					   if ( isfile(join(pathpref, f)) & \
							('.nc' in f ) & (filepref in f) ) ])

	return(filelist)

# ========================================================================================
# ========================================================================================
#                                        Main
# ========================================================================================
# ========================================================================================

def main(argv):

	# Interpret the arguments to the script
	eval_id, test_h_pref, var_name = interp_args(argv)

	# Close all figures
	plt.close('all')


        # Simple date and time
        moname  = ['Jan','Feb','Mar','Apr','May','Jun', \
                   'Jul','Aug','Sep','Oct','Nov','Dec']
        modays  = [31,28,31,30,31,30,31,31,30,31,30,31]
        cmodays = np.cumsum(modays)

	plotfile_name = eval_id+"_tseries.pdf"
	pdf = PdfPages(plotfile_name)


	# Get a list of files
	# -----------------------------------------------------------------------------------
	h0_list = getnclist(test_h_pref,'h0')



	# Load of the geographic coordinates and masks from the first file
	# -----------------------------------------------------------------------------------
	fpcoord = netcdf.netcdf_file(h0_list[0], 'r', mmap=False)
    
	# Load up the coordinate data
	latvec_in = fpcoord.variables['lat'].data;
	lonvec_in = fpcoord.variables['lon'].data;
    
	# Change coordinate system and create a re-order index array
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
	landfrac = fpcoord.variables['landfrac'].data[1:-1,sort_ids]

	# Save the ocean-ids
	ocean_ids = np.where(landfrac>1.0)
	landfrac[ocean_ids]=np.nan
        ocean_ids = np.where(landfrac<0.05)
        landfrac[ocean_ids]=np.nan
        #	SingleMapPlot(xv,yv,landfrac,ylgn_seq_cmap,[0.0,1.0],'Land Fraction',pdf)


	coldstat_raw = fpcoord.variables[var_name.strip()].data[:,1:-1,sort_ids]
	vrange = [1.0, 2.0]
	for itime in range(0,coldstat_raw.shape[0]):
                print("Day: {}".format(itime))
                imo=np.argmax(cmodays>=(itime+1))+1
                if(imo>1):
                        moday1=cmodays[imo-2]
                        day = (itime+1)-moday1
                else:
                        day = (itime+1)
		coldstat_itime = coldstat_raw[itime,:,:]
		coldstat_itime[np.where(coldstat_itime>100.0)]=np.nan
		coldstat_itime[ocean_ids]=np.nan
		ims=SingleMapPlot(itime,xv,yv,coldstat_itime,ylgn_seq_cmap,vrange,'Cold Status {0}-{1:2d}'.format(moname[imo-1],day))

        

        


def SingleMapPlot(itime,xv,yv,mapdata,color_map,vrange,map_title):

	fig = plt.figure()
        m = Basemap(projection='robin',lon_0=0,resolution='c')
	xmap,ymap = m(xv,yv)
	
	m.drawcoastlines()
	m.pcolormesh(xmap,ymap,np.ma.masked_invalid(mapdata),cmap=color_map,vmin=vrange[0],vmax=vrange[1])
	m.colorbar()
	m.drawparallels(np.arange(-90.,120.,30.))
	m.drawmeridians(np.arange(0.,360.,60.))
	plt.title(map_title)
	
	fig.savefig("frames/coldstat_frame_{}.png".format(str(itime).zfill(3)),dpi=150,frameon=False)
	plt.close(fig)








# =======================================================================================
# This is the actual call to main
   
if __name__ == "__main__":
    main(sys.argv)
