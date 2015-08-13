#!/usr/bin/env python
'''
Sifting script that knows how to take several lpps_search_script.py output
directories and sift all the 'raw' candidates. 
'''
import optparse
import os
import shutil
import random
import sys
import time

from pulpsearch import duplicates
from pulpsearch import crawler
from pulpsearch import inf
from pulpsearch import sift
from pulpsearch.util import parse_list_of_directories, parse_list_of_files
from fold import fold_candidates

#N_CANDIDATES_CUTOFF = 200
LOW_CANDIDATE = 0
HIGH_CANDIDATE = 200


def create_mixed_dirs(search_out_dirs, work_dirs):
    for wd in work_dirs:
        if len(os.listdir(wd)) > 0:
            raise Exception('Output directory %s not empty.' % wd)

    for sd, wd in zip(search_out_dirs, work_dirs):
        mixed_dir = os.path.join(wd, 'CANDIDATES')
        inf_dir = os.path.join(sd, 'INF')
        cand_dir = os.path.join(sd, 'ACCELSEARCH')

        os.makedirs(mixed_dir)
        for f in os.listdir(inf_dir):
            f = os.path.join(inf_dir, f)
            if os.path.isfile(f):
                shutil.copy(f, mixed_dir)
        for f in os.listdir(cand_dir):
            f = os.path.join(cand_dir, f)
            if os.path.isfile(f):
                shutil.copy(f, mixed_dir)
        mixed_dirs.append(mixed_dir)
    return mixed_dirs

def n_beam_sift(mixed_dirs, basenames, *args, **kwargs):
    sifted_per_beam = []
    unsifted_per_beam = []

    all_sifted_candidates = [] 
    print "basenames 2",basenames
    for md, basename in zip(mixed_dirs, basenames): 
        # assumes DM=0 to be present
        print "-----basename ",basename
        metadata = inf.inf_reader(os.path.join(md, basename + '_DM0.00.inf'))
        u, s = sift.sift_accel_cands(md, basename, metadata=metadata, n_candidates_cutoff=0)
        all_sifted_candidates.extend(s)
        unsifted_per_beam.append(u)

    n_before = len(all_sifted_candidates)
    deduplicated = duplicates.remove_duplicates(all_sifted_candidates)
    print 'Removing duplicates across beams removed %d candidates' % \
        (n_before - len(deduplicated))

    tmp = dict([(md, list()) for md in mixed_dirs])
    for c in deduplicated:
        dest = os.path.split(c.filename)[0] 
        tmp[dest].append(c)

    for md in mixed_dirs:
        sifted_per_beam.append(tmp[md])

    return unsifted_per_beam, sifted_per_beam 

if __name__ == '__main__':
    parser = optparse.OptionParser()
    parser.add_option('--o', help='Output directories.', dest='sift_out_dirs',
        metavar='OUTPUT_DIRECTIES', type='string', default='')
    parser.add_option('--i_searchout', help='LPPS style search output directories',
        metavar='SEARCH_OUT_DIRS', type='string', default='', 
        dest='search_out_dirs')
    parser.add_option('--i_mixed', help='Mixed candidate/inf directories',
        metavar='MIXED_DIRS', type='string', default='',
        dest='mixed_dirs')
    parser.add_option('--nf', help='Don\'t perform a fold.', dest='no_fold',
        metavar='NO_FOLD', action='store_true', default=False)
    parser.add_option('--ncores', metavar='N_CORES', default=8, type='int', 
        help='Number of processor cores to use default is 8.', dest='ncores')
    parser.add_option('--zap_file', type='string', metavar='ZAPFILE',
        help='Apply zaplist file (birdies) during sifting', default='',
        dest='zap_file')
    parser.add_option('--i', help='Directories containing input SUBBAND or FITS files.',
        type='string', metavar='SUBB_DIRS', dest='raw_dirs')
    parser.add_option('--masks', help='Rfifind mask files.', type='string',
        metavar='RFI_MASK_FILES', dest='masks', default='')
    parser.add_option('--filetype', help='Type of input SUBBAND (default) or FITS.',
        metavar='FILETYPE', type='string', dest='filetype')

    options, args = parser.parse_args()
    start_time = time.time()

    print 'N BEAM SIFTING / FOLDING SCRIPT, CHECKING COMMAND LINE OPTIONS.'

    # check commandline options for consistency
    if not (options.mixed_dirs or options.search_out_dirs):
        print 'Please specify either the --i_searchout or the --i_mixed option.'
        parser.print_help()
        sys.exit()

    if options.raw_dirs:
        raw_dirs = parse_list_of_directories(options.raw_dirs, 
            must_exist=True)
    else:
        raw_dirs = []

    if options.sift_out_dirs:
        sift_out_dirs = parse_list_of_directories(options.sift_out_dirs,
            must_exist=True)
    else:
        print 'No output directories provided.'
        parser.print_help()
        sys.exit()

    if options.mixed_dirs:
        mixed_dirs = parse_list_of_directories(options.mixed_dirs, 
            must_exist=True)
        search_out_dirs = []
        if len(mixed_dirs) != len(sift_out_dirs):
            print 'Specify equally many input as output directories.'
            parser.print_help()
            sys.exit()
        if raw_dirs and len(raw_dirs) != len(mixed_dirs):
            print 'Specify equally many input as output directories.'
            parser.print_help()
            sys.exit()
        basenames = []
        for md in mixed_dirs:
            basename = crawler.get_basename(md, r'^(?P<basename>\S+)_DM\d+\.\d{2}\.inf$')
            basenames.append(basename)
    else:
        search_out_dirs = parse_list_of_directories(options.search_out_dirs,
            must_exist=True)
        mixed_dirs = []
        if len(search_out_dirs) != len(sift_out_dirs):
            print 'Specify equally many input as output directories.'
            parser.print_help()
            sys.exit()
        if raw_dirs and len(raw_dirs) != len(sift_out_dirs):
            print'Specify equally many input as output directories.'
            parser.print_help()
            sys.exit()

        basenames = []
        print search_out_dirs
        for sd in search_out_dirs:
            tmp = os.path.join(sd, 'INF')
            basename = crawler.get_basename(tmp, r'^(?P<basename>\S+)_DM\d+\.\d{2}\.inf$')
            basenames.append(basename)
 
    if not options.no_fold and not raw_dirs:
        print 'If you want to fold, you need to provide the --i option.'
        parser.print_help()
        sys.exit()

    # deal with rfifind masks
    if mixed_dirs and options.masks:
        masks = parse_list_of_files(options.masks, True)
        if not options.no_fold and len(masks) != len(mixed_dirs):
            print 'Specify an rfimask for each beam.'
            parser.print_help()
            sys.exit()    
    elif search_out_dirs:
        masks = []
        for i, d in enumerate(search_out_dirs):
            files = os.listdir(os.path.join(d, 'RFIFIND'))
            L = len(basenames[i])
            found = False
            for f in files:
                if f[L:] == '_rfifind.mask':
                    masks.append(os.path.join(d, 'RFIFIND', f))
                    found = True
                    break
            if not found:
                print 'Something is wrong with search output directory %s' % d
                print 'Cannot find the rfifind mask.'
                parser.print_help()
                sys.exit(1)
 
                
    else:
        print 'Not using rfifind masks during folding, proceeding!'
        masks = []

    # deal with .zaplist files
    if options.zap_file:
        zap_file = os.path.abspath(options.zap_file)
        if not os.path.exists(zap_file):
            print 'File specified with --zap_file does not exist'
            parser.print_help()
            sys.exit(1)
    else:
        zap_file = ''

    if not options.filetype in ['SUBBAND', 'FITS', 'FIL']:
        print 'File type %s not known, exiting!\n' % options.filetype
        parser.print_help()
        sys.exit(1)

    # start processing:
    if search_out_dirs:
        print 'START COPYING FILES OVER (WILL TAKE A WHILE)'
        mixed_dirs = create_mixed_dirs(search_out_dirs, sift_out_dirs)
    print 'COMMAND LINE CHECKS OUT, PROCEEDING.'


    # perform the sift across beams
    print 'SIFTING CANDIDATES IN ALL DATA SETS (BEAMS)'
    print "basenames",basenames
    u, s = n_beam_sift(mixed_dirs, basenames, zaplist_file=zap_file)
    u_count = sum(len(x) for x in u)
    s_count = sum(len(x) for x in s)
    print 'TOTAL NUMBER OF CANDIDATES BEFORE SIFTING %d' % u_count
    print 'TOTAL NUMBER OF CANDIDATES AFTER SIFTING %d' % s_count
    
    filetype = 'SUBBAND' # <- PLACEHOLDER FOR EVENTUAL FITS FILE HANDLING TODO FIXME
    for i in range(len(raw_dirs)):
        sifted_candidates = s[i][:HIGH_CANDIDATE]
        unsifted_candidates = u[i][:HIGH_CANDIDATE]
        raw_dir = raw_dirs[i]
        work_dir = sift_out_dirs[i]
        if masks:
            mask_filename = masks[i]
        else:
            mask = '' 
        fold_candidates(sifted_candidates, raw_dir, basename, work_dir,
            options.ncores, options.filetype, mask_filename=mask_filename) 
        sift.write_candidates(os.path.join(work_dir, 'FOLDS', 'SIFTED_CANDIDATES'),
            sifted_candidates)
        # sift.write_plots(os.path.join(work_dir, 'FOLDS',  'CANIDATE_HISTOGRAMS'), 
         #   sifted_candidates, unsifted_candidates)
 
    print '\nN BEAM SIFTING / FOLDING SCRIPT IS DONE.'
    end_time = time.time()

