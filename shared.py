#!/usr/bin/env python3

import csv
import os
import subprocess
import urllib.request

def run(executable, args=[], stdin=[], stdout=None):
  pstdin = subprocess.PIPE if len(stdin) > 0 else None
  pstdout = None if stdout is None else open(stdout, "w")
  command = [executable] + args
  p = subprocess.Popen(command, stdin=pstdin, stdout=pstdout, encoding="utf8")
  if pstdin == subprocess.PIPE:
    for line in stdin:
      p.stdin.write(line + "\n")
    p.stdin.close()
  p.wait()

class Structure:
  def __init__(self, row):
    self.id = row["structureId"].lower()
    self.resolution = float(row["resolution"])
    self.rwork = float(row["rWork"])
    self.rfree = float(row["rFree"])
    self.chains = []

class Chain:
  def __init__(self, row):
    self.id = row["chainId"]
    self.cluster95 = row["clusterNumber95"]
    self.cluster90 = row["clusterNumber90"]
    self.cluster70 = row["clusterNumber70"]
    self.cluster50 = row["clusterNumber50"]
    self.cluster40 = row["clusterNumber40"]
    self.cluster30 = row["clusterNumber30"]

def get_structure_dict():
  columns = [
    "experimentalTechnique",
    "resolution",
    "rWork",
    "rFree",
    "entityMacromoleculeType",
    "clusterNumber95",
    "clusterNumber90",
    "clusterNumber70",
    "clusterNumber50",
    "clusterNumber40",
    "clusterNumber30",
  ]
  path = os.path.join("data", "pdb-chains.csv")
  if not os.path.exists(path):
    print("Downloading RCSB custom report ...")
    url = "https://www.rcsb.org/pdb/rest/customReport.xml?pdbids=*&"
    url += "customReportColumns=%s&" % ",".join(columns)
    url += "format=csv&service=wsfile"
    urllib.request.urlretrieve(url, path)
  print("Reading structures from pdb-chains.csv ...")
  structures = {}
  with open(path) as f:
    for row in csv.DictReader(f):
      if any(row[column] == "" for column in columns): continue
      if row["experimentalTechnique"] != "X-RAY DIFFRACTION": continue
      if row["entityMacromoleculeType"] != "Polypeptide(L)": continue
      structureId = row["structureId"]
      if structureId not in structures:
        structures[structureId] = Structure(row)
      structures[structureId].chains.append(Chain(row))
  return structures

def get_structure_list():
  return list(get_structure_dict().values())
