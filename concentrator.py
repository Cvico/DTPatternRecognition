''' Script to make concentrator studies '''
import os
from optparse import OptionParser
import ROOT as r
import numpy as np

from utils.ntuple_reader import ntuple
from utils.functions import color_msg
import filters
from utils.baseHistos import histos

def addConcentratorOptions(pr):
  pr.add_option('--inpath', '-i', type="string", dest = "inpath", default = "./results/")
  
  # Additional
  pr.add_option("--outfolder", "-o", type="string", dest = "outfolder", default = "./results")
  pr.add_option('--maxfiles', type=int, dest = "maxfiles", default = -1)
  pr.add_option('--maxevents', type=int, dest = "maxevents", default = -1)
  return

if __name__ == "__main__":
  pr = OptionParser(usage="%prog [options]")
  addConcentratorOptions(pr)
  (options, args) = pr.parse_args()
  inpath = options.inpath
  outfolder = options.outfolder
  maxfiles = options.maxfiles
  maxevents = options.maxevents
  
  # Analyses to be run
  # (postfix, filters)
  run_over = [
    ("_AM_withShowers", [filters.baseline]),
    ("_AM_noShowers",   [filters.baseline, filters.removeShower])
  ]
  
  for parameters in run_over:
    color_msg(f"Shower performance analyzer: {parameters[0]} ", "green")
    ntuplizer = ntuple(
      inputFolder = inpath, 
      selectors = parameters[1],
      histograms = histos,
      outfolder = outfolder, 
      maxevents = maxevents, 
      maxfiles = maxfiles,
      postfix = parameters[0]
    )
    ntuplizer.run()