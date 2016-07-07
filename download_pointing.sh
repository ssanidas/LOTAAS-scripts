#! /bin/bash

# Verbose mode
# set -x

# Read command line
IDS="$@"

# Set up Grid credentials
module load globus
echo P4ls@rgrid | voms-proxy-init -pwstdin -valid 168:00 -dont-verify-ac -voms lofar:/lofar/pulsar
GSIDIR=gsiftp://gridftp.grid.sara.nl/pnfs/grid.sara.nl/data/lofar/pulsar/tape/lotaas
SRMDIR=srm://srm.grid.sara.nl:8443/pnfs/grid.sara.nl/data/lofar/pulsar//tape/lotaas

BASEDIR=/projects/0/lotaas/data/raw/


# Loop over OBSID, SAP, BEAM
for ID in $IDS ; do 
    # Make dir
    LOCALDIR=${BASEDIR}/${ID}_red/
    mkdir -p ${LOCALDIR}
    
    # Download data
    for SAP in 0 1 2; do 
    #for SAP in 2 ; do
	for BEAM in {0..73}; do
	    
	    FILE=${ID}_SAP${SAP}_BEAM${BEAM}
#	    if [ -d /projects/lotaas/data/out/$1/SAP${SAP}/$1_SAP${SAP}_BEAM${BEAM}/BEAM${BEAM}_sift/FOLDS/CORE_0/ ]
 #           then
		#echo "dir exists"
#		if [ `ls /projects/lotaas/data/out/$1/SAP${SAP}/$1_SAP${SAP}_BEAM${BEAM}/BEAM${BEAM}_sift/FOLDS/CORE_0/*.pfd | wc -l` -gt 1 ]
#		then            
#		    echo "pfd's exist for beam ${BEAM}"
#		    continue
#		else
#		    echo "no pfd's for beam ${BEAM}"
#		fi
#	    else
#		echo "beam ${BEAM} has not been processed"
 #           fi

            #check if rfi files exist already
#	    if [ -e /projects/lotaas/data/raw/$1_red/stokes/SAP${SAP}/BEAM${BEAM}/$1_SAP${SAP}_BEAM${BEAM}_rfifind.rfi ]
 #           then
#		echo "there is a tar file"
#	    else
#		echo "no rfi file"
	        # Download and untar meta data (mask, etc).
	    globus-url-copy -vb ${GSIDIR}/${ID}/${FILE}.meta.tar file://${LOCALDIR}/${FILE}.meta.tar
	    (cd ${LOCALDIR} && tar xvf ${FILE}.meta.tar)
	    #fi
	    
	    # check if FITS exists
	    FITSDIR=${LOCALDIR}/stokes/SAP${SAP}/BEAM${BEAM}/
#	    if [ -e ${FITSDIR}${FILE}.fits ]
#	    then
#		if [ ! -s ${FITSDIR}${FILE}.fits ]
#		then
#		    echo "fits has no size - try downloading"
	    mkdir -p $FITSDIR
	    globus-url-copy -vb ${GSIDIR}/${ID}/${FILE}.fits file://${FITSDIR}/${FILE}.fits
#	        else
#		    echo "fits exists"
#		fi
#	    else
#		echo "no fits file"
#		# Download FITS file
#	    	mkdir -p $FITSDIR
#		globus-url-copy -vb ${GSIDIR}/${ID}/${FILE}.fits file://${FITSDIR}/${FILE}.fits
#	    fi
	    # Fix coordinates
	    #python /home/joerivl/bin/fix_fits_coordinates.py \
		#-v -d /projects/lotaas/parsets/ ${FITSDIR}/${FILE}.fits

#	    ssh cooper@ui.grid.sara.nl "srm-release-files ${SRMDIR}/${ID}/${FILE}.fits"
#	    ssh cooper@ui.grid.sara.nl "srm-release-files ${SRMDIR}/${ID}/${FILE}.meta.tar"
	    #ssh cooper@ui.grid.sara.nl 'srm-release-files ${SRMDIR}/${ID}/${FILE}.meta.tar'
	   
#Fix the coordinates in the coherent beams 
python /home/sanidas/scripts_lotaas/fix_fits_coordinates.py -v -d ${LOCALDIR}/ ${FITSDIR}/${FILE}.fits

#Fix the incoherent beam issue
if [[ $BEAM == 6 ]];then
    python /home/sanidas/scripts_lotaas/fix_coords_IS_beams.py --meta ${LOCALDIR}/incoherentstokes/SAP${SAP}/BEAM${BEAM}/$1_SAP*.h5 --verbose ${FITSDIR}/${FILE}.fits
fi
	done
    done
done
