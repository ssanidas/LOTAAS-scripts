#!/bin/bash

obsfile=$1
realfft -outdir . $obsfile
rm -f $obsfile
