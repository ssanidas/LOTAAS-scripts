#!/usr/bin/env python
#
#
"""
	Executable fix_coords_IS_beams.py
	For IS beams of LOTAAS pointings starting from July 1, L233346
	till L253410 (start of Dec)
	based on SAP Group metadata of HDF5 file

	(c) Vlad Kondratiev, Dec 9, 2014
"""

import os, sys, math, re
import numpy as np
import optparse as opt
import pyfits
import ephem
import h5py

# snippet from Joeri's fix_fits_coordinates.py
def write_header_comment(string, fits):
	start=fits.tell()
	fits.seek(-49,1) # go back to start
	fits.write(string)
	fits.seek(start)

# snippet from Joeri's fix_fits_coordinates.py
def write_header_pos(string, fits):
	# make sure first zero is in there
	string = re.sub("^(-?)(\d:)", "\g<1>0\g<2>", string)

	# format with quotes and some blanking spaces
	fitsval = "'%s'   " % string

	# overwrite current value
	start=fits.tell()
	fits.seek(-70,1) # go back to start
	fits.write(fitsval)
	fits.seek(start)

# read SAP coordinates from .h5 file
def get_SAP_coordinates(h5file, sapid):
	try:
		f5 = h5py.File(h5file, 'r')
		rarad = np.float64((f5['SUB_ARRAY_POINTING_%03d' % (sapid)].attrs['POINT_RA'] * math.pi) / 180.)
		decrad = np.float64((f5['SUB_ARRAY_POINTING_%03d' % (sapid)].attrs['POINT_DEC'] * math.pi) / 180.)
		f5.close()
	except:
		print "Can't get SAP%d coordinates from the %s file!" % (sapid, h5file)
		sys.exit(1)
	return (rarad, decrad)


# Main
if __name__=="__main__":
	#
	# Parsing the command line options
	#
	usage = "Usage: %prog [options] .fits"
	cmdline = opt.OptionParser(usage)
	cmdline.add_option('--meta', dest='metafile', metavar='.h5', help="HDF5 .h5 file with the meta info about the observation", default="", type='str')
	cmdline.add_option('-v', '--verbose', dest='is_verbose', action="store_true", help="Verbose output, print extra info", default=False)

	# reading cmd options
	(opts,args) = cmdline.parse_args()

	# check if input file is given
	if len(args) == 0:
		cmdline.print_usage()
		sys.exit(0)

	# reading input fits-file
	infile = args[0]
	
	# checking if meta file is given
	if opts.metafile == "":
		print "No .h5 file is given with --meta option. Can't get coordinates of the SAP!"
		sys.exit(1)

	# checking that ObsID of the fits file is within the range that should be corrected
	try:
		obsid = int(infile.split("_")[0][1:])
		if obsid < 233346 or obsid > 253408:
			print "The ObsID of the given fits-file is out the range of those that should be corrected!"
			sys.exit(1)
	except:
		print "Wrong naming of the fits-file!"
		sys.exit(1)

	# checking that given fits-file is for IS beam (BEAM12)
	try:
		beamid = int(infile.split("_BEAM")[-1].split(".fits")[0])
		if beamid != 12:
			print "This is not a fits-file for the IS beam!"
			sys.exit(1)
	except:
		print "Wrong naming of the fits-file!"
		sys.exit(1)

	# getting SAP number of the fits file
	try:
		sapid = int(infile.split("_SAP")[-1].split("_")[0])
	except:
		print "Wrong naming of the fits-file!"
		sys.exit(1)

	if (opts.is_verbose):
		print "Working on file %s:" % (infile)

	# reading SAP coordinates
	(rarad, decrad) = get_SAP_coordinates(opts.metafile, sapid)

	# convert to HHMMSS
	src = ephem.FixedBody()
	src._ra  = float(rarad)
	src._dec = float(decrad)

	Update_RA=Update_Dec=False

	fits = open(infile, "r+b")
	while not (Update_RA and Update_Dec):
		data = fits.read(80)
        
		# Fix the RA
		if (data.split()[0]=='RA'):
			if (opts.is_verbose):
				print "  Replacing RA : %s -> %s'" % (data.split()[2], str(src._ra))
			write_header_pos(str(src._ra), fits)
			write_header_comment('/ Right ascension (hh:mm:ss.ssss) [Updated]', fits)
			Update_RA=True

		# Fix the DEC
		if (data.split()[0]=='DEC'):
			if (opts.is_verbose):
				print "  Replacing DEC: %s -> %s'" % (data.split()[2], str(src._dec))
			write_header_pos(str(src._dec), fits)
			write_header_comment('/ Declination (-dd:mm:ss.sss) [Updated]', fits)
			Update_Dec=True
 
	fits.close()
