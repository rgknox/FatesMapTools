# =======================================================================================
#
# This script hold various map plotting functions that are called by the various
# other scripts in FatesMapTools.py
#
# =======================================================================================

import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import sys
import math
import code  # For development: code.interact(local=dict(globals(), **locals()))
from mpl_toolkits.basemap import Basemap

pdf_type = 1
png_type = 2

# =======================================================================================

## get a list of files given a directory prefix and specification on the type
# @param file_prefix a string with the full or relative path to the data
# @return filelist a list of strings, each of which a file

def GetNCList(file_prefix,filetype):

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


# Define the map class.  These classes can be appended into a
# list in the map_plot_type class, or they can be used stand-alone

class map_type:

	def __init__(self,data,title,co_map,var_range):

		self.data = data           # This should be a 2D rectilinear map
		self.title = title         # This is the title of the plot on this axis
		self.co_map = co_map       # This is the desired colormap
		self.var_range = var_range # This is either the desired range (hi - low)
                                   # Or, if this has more than two entries
                                   # it is the actual values to use in the
                                   # colorbar

class map_plot_title_type:

	def __init__(self,title_str,xoff,yoff,fontsize):

		self.title_str = title_str
		self.xoff      = xoff
		self.yoff      = yoff
		self.fontsize  = fontsize


# Define the plot class (holds the maps, must be the same grids)

class map_plot_type:

	# The constructor

	def __init__(self,xv,yv,proj_type,outfile_name,do_save,save_type,save_obj):

		self.xv = xv                     # The x grid coordinate
		self.yv = yv                     # The y grid coordinate
		self.proj_type = proj_type       # This is the string name of a projection
                                                 # that is found in basemap
		self.outfile_name = outfile_name # The name of the output file
		self.do_save      = do_save

                # This is either an integer equal to png_type or pdf_type
                self.save_type    = save_type

                # This is just an empty list for
                self.save_obj     = save_obj


# Call the plotting functions
def PlotMaps(plot_obj, map_list, title_obj):

	fig = plt.figure()

	# Determine the number of axes in this plot

	nplots = len(map_list)

	# If the plot title is not blank, we will add this to the

	if(len(title_obj.title_str)>0):
		npanels = nplots+1
		title_pan = True
	else:
		npanels = nplots
		title_pan = False

	if(npanels==1):
		ax_code = 110
	elif(npanels==2):
		ax_code = 120
	elif(npanels<5):
		ax_code = 220
	elif(npanels<7):
		ax_code = 320
	else:
		print('You are printing more axes than the')
		print('generic map plotting function can handle.')
		exit(2)

	# This is our axis counter, it goes up  :P
	ax_id=0


	# If we have a title panel, write it out, in spot 1
	if(title_pan):
		ax_id=ax_id+1
		ax = fig.add_subplot(ax_code+ax_id)
		txt = ax.text(title_obj.xoff,title_obj.yoff,title_obj.title_str, \
					  horizontalalignment='left',verticalalignment='center', \
					  fontsize=title_obj.fontsize)
		txt.set_clip_on(True)
		ax.axis('off')

        # Add meridion lines. If a single plot, use dx 20
        # If more than 1, use dx of 30

        if(nplots==1):
                dpar = 20.0
                dmer = 20.0
        else:
                dpar = 30.0
                dmer = 30.0


        nmer = math.floor(360.0/dmer)
        npar = math.floor(180.0/dpar)
        meridians = np.linspace(-180+0.5*(180.0-(nmer*dmer)),180-0.5*(180.0-(nmer*dmer)),nmer+1)
        parallels = np.linspace(-90+0.5*(90.0-(npar*dpar)),90-0.5*(90.0-(npar*dpar)),npar+1)


	for imap in range(0,nplots):

		# Increment the axis index
		ax_id=ax_id+1

		# Process the range of the data and how it should
		# be constrained on the map
		if(len(map_list[imap].var_range)==0):
			vrange = [np.nanmin(map_list[imap].data), \
					  np.nanmax(map_list[imap].data)]
		elif(len(map_list[imap].var_range)==1):
			print('You specified a range for a map that has only')
			print('one value. Thats weird.')
			exit(2)
		elif(len(map_list[imap].var_range)==2):
			vrange = map_list[imap].var_range
		elif(len(map_list[imap].var_range)>2):
			vrange= [map_list[imap].var_range[0], map_list[imap].var_range[-1]]


		# Generate the axis
		ax = fig.add_subplot(ax_code+ax_id)

		# Set the title
		ax.set_title(map_list[imap].title)

		# Write the map and setup the projection
                if(plot_obj.proj_type=='robin'):
                        m = Basemap(projection=plot_obj.proj_type,lon_0=0,resolution='c')
                elif(plot_obj.proj_type=='cyl'):
                        m = Basemap(projection=plot_obj.proj_type)
                else:
                        print('An unsupported projection was passed')
                        print('FatesMapFunctions.PlotMaps')
                        print('Valid options: robin,cyl ')
                        sys.exit(2)

		xmap,ymap = m(plot_obj.xv,plot_obj.yv)

		# Add some coastlines
		m.drawcoastlines()



                m.drawparallels(parallels,labels=[False,False,False,False])
                m.drawmeridians(meridians,labels=[False,False,False,False])

		# Add the data
		m.pcolormesh(xmap,ymap,np.ma.masked_invalid(map_list[imap].data), \
					 cmap=map_list[imap].co_map,vmin=vrange[0],vmax=vrange[1])

		# Add the colorbar
                if(len(map_list[imap].var_range)>2):
                        m.colorbar(ticks=map_list[imap].var_range)
                else:
                        m.colorbar()


	if(plot_obj.do_save):
                if(plot_obj.save_type == png_type):
                        fig.savefig(plot_obj.outfile_name,dpi=150,frameon=False)
                elif(plot_obj.save_type == pdf_type):
                        plot_obj.save_obj.savefig(fig)
                else:
                        print('Unknown save type')
                        exit(2)
	else:
		plt.show()


	plt.close(fig)

def DiscreteCubeHelix(N):

    base = plt.cm.get_cmap('cubehelix')
    color_list = base(np.random.randint(0,high=255,size=N))
    cmap_name = base.name + str(N)
    return base.from_list(cmap_name, color_list, N)
