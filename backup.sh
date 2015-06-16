#!/bin/bash

#This script rsyncs critical directories in a backup location

backupdir="/home/sanidas/backup"

#rsyncing the TarLog
rsync -ah /projects/0/lotaas/data/out/TarLog $backupdir