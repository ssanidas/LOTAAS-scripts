#!/usr/bin/env python
#
from re import search
from os import stat
from sys import exit
from glob import glob
from infodata import infodata
from bisect import insort
from scipy import std, stats
from scipy.signal import detrend, convolve
from ppgplot import pgopen, pgsvp, pgpap, pgswin, pgsch, pgend, pgiden, pgmtxt, pgpt, pgbin, pgbox
from numpy import multiply, int32, float32, float64, zeros, sqrt, loadtxt, fromfile, asarray, unique, nonzero, newaxis, arange, flatnonzero, ones, size, log10
from numpy import min as nummin
from numpy import max as nummax
from presto import rfft
from optparse import OptionParser

class candidate:
    def __init__(self, DM, sigma, time, bin, downfact, dt):
        self.DM = DM
        self.sigma = sigma
        self.time = time
        self.bin = bin
        self.downfact = downfact
	self.dt = dt
    def __str__(self):
        return "%7.2f %7.2f %13.6f %10d     %3d   %g\n"%\
               (self.DM, self.sigma, self.time, self.bin, self.downfact, self.dt)
    def __cmp__(self, other):
	# Sort by time (i.e. bin) by default)
        return cmp(self.bin, other.bin)

def cmp_sigma(self, other):
    #Comparison function to sort candidates by significance
    retval = -cmp(self.sigma, other.sigma)
    return retval

def fft_convolve(fftd_data, fftd_kern, lo, hi):
    """
    fft_convolve(fftd_data, fftd_kern, lo, hi):
        Perform a convolution with the complex floating point vectors
            'fftd_data' and 'fftd_kern'.  The returned vector will start at
            at bin 'lo' (must be an integer), and go up to but not
            include bin 'hi' (also an integer).
    """
    # Note:  The initial FFTs should be done like:
    # fftd_kern = rfft(kernel, -1)
    # fftd_data = rfft(data, -1)
    prod = multiply(fftd_data, fftd_kern)
    prod.real[0] = fftd_kern.real[0] * fftd_data.real[0]
    prod.imag[0] = fftd_kern.imag[0] * fftd_data.imag[0]
    return rfft(prod, 1)[lo:hi].astype(float32)

def make_fftd_kerns(downfacts, fftlen):
    fftd_kerns = []
    for downfact in downfacts:
        kern = zeros(fftlen, dtype=float32)
        # These offsets produce kernels that give results
        # equal to scipy.signal.convolve
        if downfact % 2:  # Odd number
            kern[:downfact/2+1] += 1.0
            kern[-(downfact/2):] += 1.0
        else:             # Even number
            kern[:downfact/2+1] += 1.0
            if (downfact > 2):
                kern[-(downfact/2-1):] += 1.0
        # The following normalization preserves the
        # RMS=1 characteristic of the data
        fftd_kerns.append(rfft(kern / sqrt(downfact), -1))
    return fftd_kerns

def prune_related1(hibins, hivals, downfact):
    # Remove candidates that are close to other candidates
    # but less significant.  This one works on the raw 
    # candidate arrays and uses the single downfact
    # that they were selected with.
    toremove = set()
    for ii in range(0, len(hibins)-1):
        if ii in toremove:  continue
        xbin, xsigma = hibins[ii], hivals[ii]
        for jj in range(ii+1, len(hibins)):
            ybin, ysigma = hibins[jj], hivals[jj]
            if (abs(ybin-xbin) > downfact/2):
                break
            else:
                if jj in toremove:
                    continue
                if (xsigma > ysigma):
                    toremove.add(jj)
                else:
                    toremove.add(ii)
    # Now zap them starting from the end
    toremove = sorted(toremove, reverse=True)
    for bin in toremove:
        del(hibins[bin])
        del(hivals[bin])
    return hibins, hivals
    
def prune_related2(dm_candlist, downfacts):
    # Remove candidates that are close to other candidates
    # but less significant.  This one works on the candidate 
    # instances and looks at the different downfacts of the
    # the different candidates.
    toremove = set()
    for ii in range(0, len(dm_candlist)-1):
        if ii in toremove:  continue
        xx = dm_candlist[ii]
        xbin, xsigma = xx.bin, xx.sigma
        for jj in range(ii+1, len(dm_candlist)):
            yy = dm_candlist[jj]
            ybin, ysigma = yy.bin, yy.sigma
            if (abs(ybin-xbin) > max(downfacts)/2):
                break
            else:
                if jj in toremove:
                    continue
                prox = max([xx.downfact/2, yy.downfact/2, 1])
                if (abs(ybin-xbin) <= prox):
                    if (xsigma > ysigma):
                        toremove.add(jj)
                    else:
                        toremove.add(ii)
    # Now zap them starting from the end
    toremove = sorted(toremove, reverse=True)
    for bin in toremove:
        del(dm_candlist[bin])
    return dm_candlist

def prune_border_cases(dm_candlist, offregions):
    # Ignore those that are locate in a half-width
    # of the boundary between data and padding
    #print offregions
    toremove = set()
    for ii in range(len(dm_candlist))[::-1]:
        cand = dm_candlist[ii]
        loside = cand.bin-cand.downfact/2
        hiside = cand.bin+cand.downfact/2
        if hiside < offregions[0][0]: break
        for off, on in offregions:
            if (hiside > off and loside < on):
                toremove.add(ii)
    # Now zap them starting from the end
    toremove = sorted(toremove, reverse=True)
    for ii in toremove:
        del(dm_candlist[ii])
    return dm_candlist

full_usage = """
usage:  single_pulse_search.py [options] .dat files _or_ .singlepulse files
  [-h, --help]        : Display this help
  [-m, --maxwidth]    : Set the max downsampling in sec (see below for default)
  [-p, --noplot]      : Look for pulses but do not generate a plot
  [-t, --threshold]   : Set a different threshold SNR (default=5.0)
  [-x, --xwin]        : Don't make a postscript plot, just use an X-window
  [-s, --start]       : Only plot events occuring after this time (s)
  [-e, --end]         : Only plot events occuring before this time (s)
  [-g, --glob]        : Use the files from these glob expressions (in quotes)
  [-f, --fast]        : Use a less-accurate but much faster method of detrending
  [--oneout]          : Store search results in one single .singlepulse file for all DMs
  [--dms]             : Only plot events with DM larger than this value (in pc/cc)
  [--dme]             : Only plot events with DM smaller than this value (in pc/cc)
  [--wmin]            : Only plot events with widths larger than this value (in ms)
  [--wmax]            : Only plot events with widths smaller than this value (in ms)
  [--no-width-label]  : Don't write width range label on the plot

  Perform a single-pulse search (or simply re-plot the results of a
  single-pulse search) on a set of de-dispersed time series (.dat
  files).

  The search attempts to find pulses by matched-filtering the data with a
  series of different width boxcar functions.  The possible boxcar sizes
  are [1, 2, 3, 4, 6, 9, 14, 20, 30, 45, 70, 100, 150, 250, 500, 1000]
  bins.  By default the boxcars <= 30 are used.  You can specify
  that the larger boxcars are used with the -m (or --maxwidth) option.

  The matched filtering (and accounting for all the possible 'phase'
  offsets of each boxcar) is accomplished by convolving the boxcars
  with the full resolution data.  'Duplicate' candidates from this
  process are filtered, leaving only the most significant.  The time
  series are initially smoothed using a piecewise linear fit to the
  data where each piece is 2000 data points long.

  If the input files are .singlepulse files, we won't actually perform
  a search, we'll only read in the output .singlepulse files and make
  a plot using the information they contain (along with the
  corresponding .inf files).

  Copyright Scott Ransom <sransom@nrao.edu>, 2005
  Updated version: Vlad Kondratiev, 2014
"""
usage = "usage: %prog [options] .dat files _or_ .singlepulse files"
    
def read_singlepulse_files(infiles, threshold, T_start, T_end, W_start, W_end):
    DMs = []
    candlist = []
    num_v_DMstr = {}
    if infiles[0].endswith(".singlepulse"):
        filenmbase = infiles[0][:infiles[0].rfind(".singlepulse")]
        if search("_DMs", filenmbase):
            dm=filenmbase[filenmbase.find("_DMs")+4:filenmbase.rfind("-")]
            filenmbase="%s_DM%s" % (filenmbase[:filenmbase.find("_DMs")], dm)                   
    else:
        filenmbase = infiles[0]
    info = infodata(filenmbase+".inf")

    for infile in infiles:
        if stat(infile)[6]:
            try:
                cands = loadtxt(infile)
                if len(cands.shape)==1:
                    cands = asarray([cands])
                for cand in cands:
                    if cand[2] < T_start: continue
                    if cand[2] > T_end: continue
                    width = cand[4] * cand[5] * 1000.
                    if width < W_start: continue
                    if width > W_end: continue
                    if cand[1] >= threshold:
                        candlist.append(candidate(*cand))
			DMs.append(float(cand[0]))
			DMstr="%.2f" % (cand[0])
			if DMstr in num_v_DMstr:
                            num_v_DMstr[DMstr] += 1
			else:
			    num_v_DMstr[DMstr] = 0
            except:  # No candidates in the file
                IndexError

    DMs=list(unique(asarray(DMs)))
#    DMs.sort()   # no need here as it will be sorted anyway before plotting
    return info, DMs, candlist, num_v_DMstr

def main():
    parser = OptionParser(usage)
    parser.add_option("-x", "--xwin", action="store_true", dest="xwin",
                      default=False, help="Don't make a postscript plot, just use an X-window")
    parser.add_option("-p", "--noplot", action="store_false", dest="makeplot",
                      default=True, help="Look for pulses but do not generate a plot")
    parser.add_option("-m", "--maxwidth", type="float", dest="maxwidth", default=0.0,
                      help="Set the max downsampling in sec (see below for default)")
    parser.add_option("-t", "--threshold", type="float", dest="threshold", default=5.0,
                      help="Set a different threshold SNR (default=5.0)")
    parser.add_option("-s", "--start", type="float", dest="T_start", default=0.0,
                      help="Only plot events occuring after this time (s)")
    parser.add_option("-e", "--end", type="float", dest="T_end", default=1e9,
                      help="Only plot events occuring before this time (s)")
    parser.add_option("-g", "--glob", type="string", dest="globexp", default=None,
                      help="Process the files from this glob expression")
    parser.add_option("-f", "--fast", action="store_true", dest="fast",
                      default=False, help="Use a faster method of de-trending (2x speedup)")
    parser.add_option("--oneout", action="store_true", dest="oneout",
                      default=False, help="Store search results in one single .singlepulse file for all DMs")
    parser.add_option("--dms", type="float", dest="DM_start", default=0.0,
                      help="Only plot events with DM larger than this value (in pc/cc)")
    parser.add_option("--dme", type="float", dest="DM_end", default=1e9,
                      help="Only plot events with DM smaller than this value (in pc/cc)")
    parser.add_option("--wmin", type="float", dest="W_start", default=0.0,
                      help="Only plot events with widths larger than this value (in ms)")
    parser.add_option("--wmax", type="float", dest="W_end", default=1000000.0,
                      help="Only plot events with widths smaller than this value (in ms)")
    parser.add_option("--no-width-label", action="store_true", dest="is_no_width_label",
                      default=False, help="Don't write width range label on the plot")
    (opts, args) = parser.parse_args()
    if len(args)==0:
        if opts.globexp==None:
            print full_usage
            exit(0)
        else:
            args = []
            for globexp in opts.globexp.split():
                args += glob.glob(globexp)
    useffts = True
    dosearch = True
    if opts.xwin:
        pgplot_device = "/XWIN"
    else:
        pgplot_device = ""

    fftlen = 8192     # Should be a power-of-two for best speed
    chunklen = 8000   # Must be at least max_downfact less than fftlen
    detrendlen = 1000 # length of a linear piecewise chunk of data for detrending
    blocks_per_chunk = chunklen / detrendlen
    overlap = (fftlen - chunklen)/2
    worklen = chunklen + 2*overlap  # currently it is fftlen...

    max_downfact = 30
    default_downfacts = [2, 3, 4, 6, 9, 14, 20, 30, 45, 70, 100, 150, 250, 500, 1000]

    if args[0].endswith(".singlepulse"):
        filenmbase = args[0][:args[0].rfind(".singlepulse")]
        dosearch = False
    elif args[0].endswith(".dat"):
        filenmbase = args[0][:args[0].rfind(".dat")]
    else:
        filenmbase = args[0]

    # Don't do a search, just read results and plot
    if not dosearch:
        info, DMs, candlist, num_v_DMstr = \
              read_singlepulse_files(args, opts.threshold, opts.T_start, opts.T_end, opts.W_start, opts.W_end)
        orig_N, orig_dt = int(info.N), info.dt
        obstime = orig_N * orig_dt
    else:
        DMs = []
        candlist = []
        num_v_DMstr = {}

        # Loop over the input files
        for filenm in args:
            if filenm.endswith(".dat"):
                filenmbase = filenm[:filenm.rfind(".dat")]
            else:
                filenmbase = filenm

            info = infodata(filenmbase+".inf")
            DMstr = "%.2f"%info.DM
            DMs.append(info.DM)
            N, dt = int(info.N), info.dt
            obstime = N * dt
            # Choose the maximum width to search based on time instead
            # of bins.  This helps prevent increased S/N when the downsampling
            # changes as the DM gets larger.
            if opts.maxwidth > 0.0:
                downfacts = [x for x in default_downfacts if x*dt <= opts.maxwidth]
            else:
                downfacts = [x for x in default_downfacts if x <= max_downfact]
            if len(downfacts) == 0:
                downfacts = [default_downfacts[0]]
            if (filenm == args[0]):
                orig_N = N
                orig_dt = dt
                if useffts:
                    fftd_kerns = make_fftd_kerns(downfacts, fftlen)
            if info.breaks:
                offregions = zip([x[1] for x in info.onoff[:-1]],
                                 [x[0] for x in info.onoff[1:]])

            # Compute the file length in detrendlens
            roundN = N/detrendlen * detrendlen
            numchunks = roundN / chunklen
            # Read in the file
            print 'Reading "%s"...'%filenm
            timeseries = fromfile(filenm, dtype=float32, count=roundN)
            # Split the timeseries into chunks for detrending
            numblocks = roundN/detrendlen
            timeseries.shape = (numblocks, detrendlen)
            stds = zeros(numblocks, dtype=float64)
            # de-trend the data one chunk at a time
            print '  De-trending the data and computing statistics...'
            for ii, chunk in enumerate(timeseries):
                if opts.fast:  # use median removal instead of detrending (2x speedup)
                    tmpchunk = chunk.copy()
                    tmpchunk.sort()
                    med = tmpchunk[detrendlen/2]
                    chunk -= med
                    tmpchunk -= med
                else:
                    # The detrend calls are the most expensive in the program
                    timeseries[ii] = detrend(chunk, type='linear')
                    tmpchunk = timeseries[ii].copy()
                    tmpchunk.sort()
                # The following gets rid of (hopefully) most of the 
                # outlying values (i.e. power dropouts and single pulses)
                # If you throw out 5% (2.5% at bottom and 2.5% at top)
                # of random gaussian deviates, the measured stdev is ~0.871
                # of the true stdev.  Thus the 1.0/0.871=1.148 correction below.
                # The following is roughly .std() since we already removed the median
                stds[ii] = sqrt((tmpchunk[detrendlen/40:-detrendlen/40]**2.0).sum() /
                                    (0.95*detrendlen))
            stds *= 1.148
            # sort the standard deviations and separate those with
            # very low or very high values
            sort_stds = stds.copy()
            sort_stds.sort()
            # identify the differences with the larges values (this
            # will split off the chunks with very low and very high stds
            locut = (sort_stds[1:numblocks/2+1] -
                     sort_stds[:numblocks/2]).argmax() + 1
            hicut = (sort_stds[numblocks/2+1:] -
                     sort_stds[numblocks/2:-1]).argmax() + numblocks/2 - 2
            std_stds = std(sort_stds[locut:hicut])
            median_stds = sort_stds[(locut+hicut)/2]
            lo_std = median_stds - 4.0 * std_stds
            hi_std = median_stds + 4.0 * std_stds
            # Determine a list of "bad" chunks.  We will not search these.
            bad_blocks = nonzero((stds < lo_std) | (stds > hi_std))[0]
            print "    pseudo-median block standard deviation = %.2f" % (median_stds)
            print "    identified %d bad blocks out of %d (i.e. %.2f%%)" % \
                  (len(bad_blocks), len(stds),
                   100.0*float(len(bad_blocks))/float(len(stds)))
            stds[bad_blocks] = median_stds
            print "  Now searching..."

            # Now normalize all of the data and reshape it to 1-D
            timeseries /= stds[:,newaxis]
            timeseries.shape = (roundN,)
            # And set the data in the bad blocks to zeros
            # Even though we don't search these parts, it is important
            # because of the overlaps for the convolutions
            for bad_block in bad_blocks:
                loind, hiind = bad_block*detrendlen, (bad_block+1)*detrendlen
                timeseries[loind:hiind] = 0.0
            # Convert to a set for faster lookups below
            bad_blocks = set(bad_blocks)

            # Step through the data
            dm_candlist = []
            for chunknum in range(numchunks):
                loind = chunknum*chunklen-overlap
                hiind = (chunknum+1)*chunklen+overlap
                # Take care of beginning and end of file overlap issues
                if (chunknum==0): # Beginning of file
                    chunk = zeros(worklen, dtype=float32)
                    chunk[overlap:] = timeseries[loind+overlap:hiind]
                elif (chunknum==numchunks-1): # end of the timeseries
                    chunk = zeros(worklen, dtype=float32)
                    chunk[:-overlap] = timeseries[loind:hiind-overlap]
                else:
                    chunk = timeseries[loind:hiind]

                # Make a set with the current block numbers
                lowblock = blocks_per_chunk * chunknum
                currentblocks = set(arange(blocks_per_chunk) + lowblock)
                localgoodblocks = asarray(list(currentblocks -
                                                   bad_blocks)) - lowblock
                # Search this chunk if it is not all bad
                if len(localgoodblocks):
                    # This is the good part of the data (end effects removed)
                    goodchunk = chunk[overlap:-overlap]

                    # need to pass blocks/chunklen, localgoodblocks
                    # dm_candlist, dt, opts.threshold to cython routine

                    # Search non-downsampled data first
                    # NOTE:  these nonzero() calls are some of the most
                    #        expensive calls in the program.  Best bet would 
                    #        probably be to simply iterate over the goodchunk
                    #        in C and append to the candlist there.
                    hibins = flatnonzero(goodchunk>opts.threshold)
                    hivals = goodchunk[hibins]
                    hibins += chunknum * chunklen
                    hiblocks = hibins/detrendlen
                    # Add the candidates (which are sorted by bin)
                    for bin, val, block in zip(hibins, hivals, hiblocks):
                        if block not in bad_blocks:
                            time = bin * dt
                            dm_candlist.append(candidate(info.DM, val, time, bin, 1, dt))

                    # Prepare our data for the convolution
                    if useffts: fftd_chunk = rfft(chunk, -1)

                    # Now do the downsampling...
                    for ii, downfact in enumerate(downfacts):
                        if useffts: 
                            # Note:  FFT convolution is faster for _all_ downfacts, even 2
                            goodchunk = fft_convolve(fftd_chunk, fftd_kerns[ii],
                                                     overlap, -overlap)
                        else:
                            # The normalization of this kernel keeps the post-smoothing RMS = 1
                            kernel = ones(downfact, dtype=float32) / \
                                     sqrt(downfact)
                            smoothed_chunk = convolve(chunk, kernel, 1)
                            goodchunk = smoothed_chunk[overlap:-overlap]
                        #hibins = nonzero(goodchunk>opts.threshold)[0]
                        hibins = flatnonzero(goodchunk>opts.threshold)
                        hivals = goodchunk[hibins]
                        hibins += chunknum * chunklen
                        hiblocks = hibins/detrendlen
                        hibins = hibins.tolist()
                        hivals = hivals.tolist()
                        # Now walk through the new candidates and remove those
                        # that are not the highest but are within downfact/2
                        # bins of a higher signal pulse
                        hibins, hivals = prune_related1(hibins, hivals, downfact)
                        # Insert the new candidates into the candlist, but
                        # keep it sorted...
                        for bin, val, block in zip(hibins, hivals, hiblocks):
                            if block not in bad_blocks:
                                time = bin * dt
                                insort(dm_candlist, candidate(info.DM, val, time, bin, downfact, dt))

            # Now walk through the dm_candlist and remove the ones that
            # are within the downsample proximity of a higher
            # signal-to-noise pulse
            dm_candlist = prune_related2(dm_candlist, downfacts)
            print "  Found %d pulse candidates"%len(dm_candlist)
            
            # Get rid of those near padding regions
            if info.breaks: prune_border_cases(dm_candlist, offregions)

            # Write the pulses to an ASCII output file
            if not opts.oneout:
                if len(dm_candlist):
                    #dm_candlist.sort(cmp_sigma)
                    outline="# DM      Sigma      Time (s)     Sample    Downfact   Sampling (s)\n"
                    for cand in dm_candlist: outline+=str(cand)
                    outfile = open(filenmbase+'.singlepulse', mode='w')
                    outfile.write(outline)
                    outfile.close()

            # Add these candidates to the overall candidate list
            for cand in dm_candlist:
		if cand.time < opts.T_start: continue
		if cand.time > opts.T_end: continue
		width = cand.downfact * cand.dt * 1000.
		if width < opts.W_start: continue
		if width > opts.W_end: continue
                candlist.append(cand)
            num_v_DMstr[DMstr] = len(candlist)

        # Writing single output .singlepulse file for a range of DMs
	if opts.oneout:
            disps=asarray([cand.DM for cand in candlist])
	    if size(disps):
	        if search("_DM", filenmbase):
                    outfilebase = filenmbase[:filenmbase.find("_DM")]
	        else:
                    outfilebase = filenmbase
                #candlist.sort(cmp_sigma)
                outline="# DM      Sigma      Time (s)     Sample    Downfact   Sampling (s)\n"
                for cand in candlist: outline+=str(cand)
	        outfile=open("%s_DMs%.2f-%.2f.singlepulse" % (outfilebase, nummin(disps), nummax(disps)), mode='w')
                outfile.write(outline)
	        outfile.close()
	

    if (opts.makeplot):

	if len(DMs):
            mindm = nummin(DMs)
            maxdm = nummax(DMs)
            DMs = [dm for dm in DMs if dm >= opts.DM_start and dm <= opts.DM_end]
            DMs.sort()
        else:
	    exit(0)  # nothing to plot

	widths = [cand.downfact*cand.dt*1000. for cand in candlist if cand.DM >= opts.DM_start and cand.DM <= opts.DM_end]
        if len(widths):
	    minw=nummin(widths)
	    maxw=nummax(widths)
        else:
	    exit(0)  # nothing to plot

        # Step through the candidates to make a SNR list
        snrs = []
        for cand in candlist:
	    if cand.DM < opts.DM_start: continue
	    if cand.DM > opts.DM_end: continue
            snrs.append(cand.sigma)
        if snrs:
            maxsnr = max(int(max(snrs)), int(opts.threshold)) + 3
        else:
            maxsnr = int(opts.threshold) + 3

        # Generate the SNR histogram
        snrs = asarray(snrs)
        (num_v_snr, lo_snr, d_snr, num_out_of_range) = \
                    stats.histogram(snrs,
                                          int(maxsnr-opts.threshold+1),
                                          [opts.threshold, maxsnr])
        snrs = arange(maxsnr-opts.threshold+1, dtype=float64) * d_snr \
               + lo_snr + 0.5*d_snr
        num_v_snr = num_v_snr.astype(float32)
        num_v_snr[num_v_snr==0.0] = 0.001

        # Generate the DM histogram
        if len(DMs):
            num_v_DM = zeros(len(DMs))
            for ii, DM in enumerate(DMs):
                num_v_DM[ii] = num_v_DMstr["%.2f"%DM]
            DMs = asarray(DMs)
        else:
            num_v_DM=[]

        # open the plot device
	short_filenmbase = filenmbase.split(".singlepulse")[0]
        if not search("_DMs", short_filenmbase):
            short_filenmbase = short_filenmbase.split("_DM")[0]

	if opts.DM_end > maxdm:
	    opts.DM_end = maxdm
        if opts.DM_start < mindm:
	    opts.DM_start = mindm

        if opts.T_end > obstime:
            opts.T_end = obstime

        if pgplot_device:
            pgopen(pgplot_device)
        else:
            if (opts.T_start > 0.0 or opts.T_end < obstime):
                if (opts.DM_start > mindm or opts.DM_end < maxdm):
                    short_filenmbase = short_filenmbase.split("_DM")[0]
                    if (opts.W_start > 0.0 or opts.W_end < 1000000.0):
                        if opts.W_end > maxw: opts.W_end = maxw
                        if opts.W_start < minw: opts.W_start = minw
                        pgopen(short_filenmbase+'_%.0f-%.0fs_DMs%.2f-%.2f_WTH%.2f-%.2fms_singlepulse.ps/VPS'%(opts.T_start, opts.T_end, opts.DM_start, opts.DM_end, opts.W_start, opts.W_end))
                    else:
                        pgopen(short_filenmbase+'_%.0f-%.0fs_DMs%.2f-%.2f_singlepulse.ps/VPS'%(opts.T_start, opts.T_end, opts.DM_start, opts.DM_end))
                else:
                    if (opts.W_start > 0.0 or opts.W_end < 1000000.0):
                        if opts.W_end > maxw: opts.W_end = maxw
                        if opts.W_start < minw: opts.W_start = minw
                        pgopen(short_filenmbase+'_%.0f-%.0fs_WTH%.2f-%.2fms_singlepulse.ps/VPS'%(opts.T_start, opts.T_end, opts.W_start, opts.W_end))
                    else:
                        pgopen(short_filenmbase+'_%.0f-%.0fs_singlepulse.ps/VPS'%(opts.T_start, opts.T_end))
            else:
                if (opts.DM_start > mindm or opts.DM_end < maxdm):
                    short_filenmbase = short_filenmbase.split("_DM")[0]
                    if (opts.W_start > 0.0 or opts.W_end < 1000000.0):
                        if opts.W_end > maxw: opts.W_end = maxw
                        if opts.W_start < minw: opts.W_start = minw
                        pgopen(short_filenmbase+'_DMs%.2f-%.2f_WTH%.2f-%.2fms_singlepulse.ps/VPS'%(opts.DM_start, opts.DM_end, opts.W_start, opts.W_end))
                    else:
                        pgopen(short_filenmbase+'_DMs%.2f-%.2f_singlepulse.ps/VPS'%(opts.DM_start, opts.DM_end))
                else:
                    if (opts.W_start > 0.0 or opts.W_end < 1000000.0):
                        if opts.W_end > maxw: opts.W_end = maxw
                        if opts.W_start < minw: opts.W_start = minw
                        pgopen(short_filenmbase+'_WTH%.2f-%.2fms_singlepulse.ps/VPS'%(opts.W_start, opts.W_end))
                    else:
                        pgopen(short_filenmbase+'_singlepulse.ps/VPS')

	if opts.W_end > maxw:
	    opts.W_end = maxw
        if opts.W_start < minw:
	    opts.W_start = minw

        pgpap(7.5, 1.0)  # Width in inches, aspect

        # plot the SNR histogram
        pgsvp(0.06, 0.31, 0.6, 0.87)
        pgswin(opts.threshold, maxsnr,
                       log10(0.5), log10(2*max(num_v_snr)))
        pgsch(0.8)
        pgbox("BCNST", 0, 0, "BCLNST", 0, 0)
        pgmtxt('B', 2.5, 0.5, 0.5, "Signal-to-Noise")
        pgmtxt('L', 1.8, 0.5, 0.5, "Number of Pulses")
        pgsch(1.0)
        pgbin(snrs, log10(num_v_snr), 1)

        # plot the DM histogram
        pgsvp(0.39, 0.64, 0.6, 0.87)
        # Add [1] to num_v_DM in YMAX below so that YMIN != YMAX when max(num_v_DM)==0
        pgswin(opts.DM_start-0.5, opts.DM_end+0.5, 0.0, 1.1*max(num_v_DM+[1]))
        pgsch(0.8)
        pgbox("BCNST", 0, 0, "BCNST", 0, 0)
        pgmtxt('B', 2.5, 0.5, 0.5, "DM (pc cm\u-3\d)")
        pgmtxt('L', 1.8, 0.5, 0.5, "Number of Pulses")
        pgsch(1.0)
        if len(DMs): pgbin(DMs, num_v_DM, 1)

        # plot the SNR vs DM plot 
        pgsvp(0.72, 0.97, 0.6, 0.87)
        pgswin(opts.DM_start-0.5, opts.DM_end+0.5, opts.threshold, maxsnr)
        pgsch(0.8)
        pgbox("BCNST", 0, 0, "BCNST", 0, 0)
        pgmtxt('B', 2.5, 0.5, 0.5, "DM (pc cm\u-3\d)")
        pgmtxt('L', 1.8, 0.5, 0.5, "Signal-to-Noise")
        pgsch(1.0)
        cand_ts = zeros(len(candlist), dtype=float32)
        cand_SNRs = zeros(len(candlist), dtype=float32)
        cand_DMs = zeros(len(candlist), dtype=float32)
        for ii, cand in enumerate(candlist):
            cand_ts[ii], cand_SNRs[ii], cand_DMs[ii] = \
                         cand.time, cand.sigma, cand.DM
        pgpt(cand_DMs, cand_SNRs, 20)

        # plot the DM vs Time plot
        pgsvp(0.06, 0.97, 0.08, 0.52)
        pgswin(opts.T_start, opts.T_end, opts.DM_start-0.5, opts.DM_end+0.5)
        pgsch(0.8)
        pgbox("BCNST", 0, 0, "BCNST", 0, 0)
        pgmtxt('B', 2.5, 0.5, 0.5, "Time (s)")
        pgmtxt('L', 1.8, 0.5, 0.5, "DM (pc cm\u-3\d)")
        # Circles are symbols 20-26 in increasing order
        snr_range = 12.0
        cand_symbols = (cand_SNRs-opts.threshold)/snr_range * 6.0 + 20.5
        cand_symbols = cand_symbols.astype(int32)
        cand_symbols[cand_symbols>26] = 26
        for ii in [26, 25, 24, 23, 22, 21, 20]:
            inds = nonzero(cand_symbols==ii)[0]
            pgpt(cand_ts[inds], cand_DMs[inds], ii)

        # Now fill the infomation area
        pgsvp(0.05, 0.95, 0.87, 0.97)
        pgsch(1.0)
        pgmtxt('T', 0.5, 0.0, 0.0,
                       "Single pulse results for '%s'"%short_filenmbase)
	if not opts.is_no_width_label:
            if opts.W_end != 1e9:
                pgmtxt('T', 0.5, 0.70, 0.0, "WIDTH=%.2f-%.2f ms" % (opts.W_start, opts.W_end))
	    else:
                pgmtxt('T', 0.5, 0.70, 0.0, "WIDTH=%.2f+ ms" % (opts.W_start))
        pgsch(0.8)
        # first row
        pgmtxt('T', -1.1, 0.02, 0.0, 'Source: %s'%\
                       info.object)
        pgmtxt('T', -1.1, 0.37, 0.0, 'RA (J2000):')
        pgmtxt('T', -1.1, 0.54, 0.0, info.RA)
        pgmtxt('T', -1.1, 0.75, 0.0, 'N samples: %.0f'%orig_N)
        # second row
        pgmtxt('T', -2.4, 0.02, 0.0, 'Telescope: %s'%\
                       info.telescope)
        pgmtxt('T', -2.4, 0.37, 0.0, 'DEC (J2000):')
        pgmtxt('T', -2.4, 0.54, 0.0, info.DEC)
        pgmtxt('T', -2.4, 0.75, 0.0, 'Sampling: %.2f \gms'%\
                       (orig_dt*1e6))
        # third row
        if info.instrument.find("pigot") >= 0:
            instrument = "Spigot"
        else:
            instrument = info.instrument
        pgmtxt('T', -3.7, 0.02, 0.0, 'Instrument: %s'%instrument)
        if (info.bary):
            pgmtxt('T', -3.7, 0.37, 0.0, 'MJD\dbary\u: %.12f'%info.epoch)
        else:
            pgmtxt('T', -3.7, 0.37, 0.0, 'MJD\dtopo\u: %.12f'%info.epoch)
        pgmtxt('T', -3.7, 0.75, 0.0, 'Freq\dctr\u: %.1f MHz'%\
                       ((info.numchan/2-0.5)*info.chan_width+info.lofreq))
        pgiden()
        pgend()

if __name__ == '__main__':
    if (0):
        # The following is for profiling
        from hotshot import Profile
        prof = hotshot.Profile("hotshot_edi_stats")
        prof.runcall(main)
        prof.close()
        # To see the results:
        if (0):
            from hotshot import stats
            s = stats.load("hotshot_edi_stats")
            s.sort_stats("time").print_stats()
    else:
        main()
