# =======================================================================================
##
# For usage: $python FatesGridPhenSeries.py
#
# This script is intended to diagnose the output of gridded CLM/ELM/FATES
# simulations, for time-series visualization of phenology variables.
#
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
from FatesMapFunctions import map_type, map_plot_type, map_plot_title_type, PlotMaps
from scipy.io import netcdf
from mpl_toolkits.basemap import Basemap




# Some constants
g_to_Mg = 1.0e-6
m2_to_ha = 1.0e-4

ylgn_seq_cmap=mpl.cm.get_cmap('YlGn')
rdbu_div_cmap=mpl.cm.get_cmap('RdBu')
wint_seq_cmap=mpl.cm.get_cmap('cool')

def usage():
	print('')
	print('=======================================================================')
	print('')
	print(' python acre_gridcomp.py -h --test-hist-pref=<text>')
	print('                         --eval-id=<text>')
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
	print(' --save-pref=<path>')
	print('     Optional.  If this path exists, instead of generating the plot')
	print('     in the window, it will save to file. The file name will have this')
	print('     prefix, with the eval-id test string appended, and then a counter')
	print('     appended to that')
	print('')
	print(' --test-hist-pref=<path>')
	print('     the full path to the history file folder')
	print('     version of output')
	print('')
	print('=======================================================================')


# ========================================================================================

## interp_args processes the arguments passed to the python script at executions.

def interp_args(argv):

	argv.pop(0)  # The script itself is the first argument, forget it

	## history file from the test simulation
	test_h_pref = ''
	## Name of the evaluation being performed, this is non-optional
	eval_id = ''
	## Name of the save location prefix (default to nothing/off)
	save_pref = ''

	try:
		opts, args = getopt.getopt(argv, 'h',["help", \
                                              "eval-id=", \
											  "save-pref=", \
                                              "test-hist-pref="])

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
		elif o in ("--save-pref"):
			save_pref = a
		elif o in ("--test-hist-pref"):
			test_h_pref = a
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


	# Remove Whitespace in eval_id string
	eval_id.replace(" ","")
    
	return (eval_id, test_h_pref, save_pref )


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
	eval_id, test_h_pref, save_pref = interp_args(argv)

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

	fpcoord.close()

	
	if(len(save_pref)>0):
		do_save = True
	else:
		do_save = False


	#code.interact(local=dict(globals(), **locals())) 
	atime=0
	for ifile in range(0,len(h0_list)):

		fpdata = netcdf.netcdf_file(h0_list[ifile], 'r', mmap=False)
				
		# Load date information
		mcdate_raw = fpdata.variables['mcdate'].data

		# Load the cold status flags, GDD and number of cold days
		coldstat_raw = fpdata.variables['SITE_COLD_STATUS'].data[:,1:-1,sort_ids]
		coldstat_raw[np.where(coldstat_raw>100.0)] = np.nan
		gdd_raw = fpdata.variables['SITE_GDD'].data[:,1:-1,sort_ids]
		gdd_raw[np.where(gdd_raw>1.e6)] = np.nan
		ncolddays_raw = fpdata.variables['SITE_NCOLDDAYS'].data[:,1:-1,sort_ids]
		ncolddays_raw[np.where(ncolddays_raw>1.e6)] = np.nan
	
		# Load the drought flags and mean liquid volume used to set status

		dstat_raw = fpdata.variables['SITE_DROUGHT_STATUS'].data[:,1:-1,sort_ids]
		dstat_raw[np.where(dstat_raw>100.0)] = np.nan

		meanliqvol_raw = fpdata.variables['SITE_MEANLIQVOL_DROUGHTPHEN'].data[:,1:-1,sort_ids]
		meanliqvol_raw[np.where(meanliqvol_raw>100.0)] = np.nan

		# Load the days since each drought and cold flags were tripped
#		cleafoff_raw = fpdata.variables['SITE_DAYSINCE_COLDLEAFOFF'].data[:,1:-1,sort_ids]
#		cleafoff_raw[np.where(cleafoff_raw>1.e6)] = np.nan
#		cleafon_raw = fpdata.variables['SITE_DAYSINCE_COLDLEAFON'].data[:,1:-1,sort_ids]
#		cleafon_raw[np.where(cleafon_raw>1.e6)] = np.nan
		
#		dleafoff_raw = fpdata.variables['SITE_DAYSINCE_DROUGHTLEAFOFF'].data[:,1:-1,sort_ids]
#		dleafoff_raw[np.where(dleafoff_raw>1.e6)] = np.nan
#		dleafon_raw = fpdata.variables['SITE_DAYSINCE_DROUGHTLEAFON'].data[:,1:-1,sort_ids]
#		dleafon_raw[np.where(dleafon_raw>1.e6)] = np.nan

		firstday=0

		for itime in range(firstday,coldstat_raw.shape[0]):

			# Absolute counter for file names
			atime=atime+1

			# Save the current date integer to a float
			# YYYYMMDD
			fdate=float(mcdate_raw[itime])

			yr   = np.floor(fdate/10000.0)
			mo   = np.floor((fdate-yr*10000.0)/100.0)
			dom  = fdate - (yr*10000.0 + mo*100.0)

			print('{} {}-{}'.format(int(yr),moname[int(mo)-1],int(dom)))

			outfile_name = "{}{}_{}.png".format(save_pref,eval_id,str(atime).zfill(4))
			proj_type    = 'robin'
			plot_title   = 'Phenology Diagnostics\n{} {}-{}'.format(int(yr),moname[int(mo)-1],int(dom))

			title_obj = map_plot_title_type(plot_title,0.25,0.4,14)
			plot_obj  = map_plot_type(xv,yv,proj_type,outfile_name,do_save)

			map_list = []

			map1=map_type(coldstat_raw[itime,:,:], 'Cold Status', wint_seq_cmap, [0,1,2])
			map_list.append(map1)
			
			map2=map_type(gdd_raw[itime,:,:], 'Growing Degree Days', ylgn_seq_cmap, [0.,100.])
			map_list.append(map2)

			map3=map_type(ncolddays_raw[itime,:,:], 'Number of Cold Days', ylgn_seq_cmap, [])
			map_list.append(map3)

			map4=map_type(dstat_raw[itime,:,:], 'Drought Status', wint_seq_cmap, [0,1,2,3])
			map_list.append(map4)
			
			map5=map_type(meanliqvol_raw[itime,:,:], 'Soil Water [m3/m3]', wint_seq_cmap, [])
			map_list.append(map5)
			


			PlotMaps(plot_obj,map_list,title_obj)

			
			

#				map3 = ncolddays_raw[itime,:,:]
#				map3[ocean_ids]=np.nan
#				map3_title = 'Number of Cold Days'
#				map3_cmap  = ylgn_seq_cmap
#				map3_vrange = [0.0,np.nanmax(map3)]

#				map4 = leafoff_raw[itime,:,:]
#				map4[ocean_ids]=np.nan
#				map4_title = 'Days Since Leaf Off'
#				map4_cmap  = ylgn_seq_cmap
#				map4_vrange = [0.0,min(400.0,np.nanmax(map4))]

#				map5 = leafon_raw[itime,:,:]
#				map5[ocean_ids]=np.nan
#				map5_title = 'Days Since Leaf On'
#				map5_cmap  = ylgn_seq_cmap
#				map5_vrange = [0.0,np.nanmax(map5)]

#				FiveMapPlot("{}{}".format('frames/',eval_id),atime,yr,moname[int(mo)-1],dom,xv,yv, \
#							map1,map2,map3,map4,map5, \
#							map1_title,map2_title,map3_title,map4_title,map5_title, \
#							map1_vrange,map2_vrange,map3_vrange,map4_vrange,map5_vrange, \
#							map1_cmap,map2_cmap,map3_cmap,map4_cmap,map5_cmap)
#			else:
#				map1 = dstat_raw[itime,:,:]
#				map1[ocean_ids]=np.nan
#				map1_title = 'Drought Status'
#				map1_cmap  = wint_seq_cmap
#				map1_vrange = [0,1,2,3]
#				SingleMapSave("{}{}".format('frames/dstatus',eval_id),atime,yr,moname[int(mo)-1],dom,xv,yv, \
#							  map1, map1_title, map1_vrange, map1_cmap)
                
#			fpdata.close()


   




# =======================================================================================
# This is the actual call to main
   
if __name__ == "__main__":
    main(sys.argv)
