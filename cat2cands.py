import logging
import subprocess
try:
  import json
  import gspread
  from oauth2client.client import SignedJwtAssertionCredentials
except ImportError:
  logging.warning("Results cannot be uploaded on the website - Module missing")
import os
import time
import argparse

from Paths import *
import Parameters


def upload_sheet(newpng):
  try: json_key = json.load(open(SITE_CERT))
  except IOError:
    logging.warning("Spreadsheet cannot be uploaded - Google certificate missing")
    return
  scope = ['https://spreadsheets.google.com/feeds']
  credentials = SignedJwtAssertionCredentials(json_key['client_email'], json_key['private_key'], scope)
  gc = gspread.authorize(credentials)
  wks = gc.open("LOTAAS_ambiguous_candidates").sheet1
  
  link = '=HYPERLINK("http://www.astron.nl/lofarpwg/lotaas/promising-candidates/{}","Plot")'.format(newpng)
  
  nrow = wks.row_count+1
  tick = wks.acell('M2').value
  countvotes = '=COUNTIF(B{}:G{}, "{}")/6'.format(nrow,nrow,tick)

  row = [link,'','','','','','',countvotes]
  wks.append_row(row)
  
if __name__=="__main__":
  parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('-f', nargs=1)
  args = parser.parse_args()

  upload_sheet(args.f[0])

