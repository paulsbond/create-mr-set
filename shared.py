#!/usr/bin/env python3

import csv
import urllib.request
import os

class Structure:
  def __init__(self, row):
    self.id = row["structureId"].lower()
    self.resolution = float(row["resolution"])
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

def get_structures():
  columns = [
    "experimentalTechnique",
    "resolution",
    "entityMacromoleculeType",
    "clusterNumber95",
    "clusterNumber90",
    "clusterNumber70",
    "clusterNumber50",
    "clusterNumber40",
    "clusterNumber30",
  ]
  if not os.path.exists("pdb-chains.csv"):
    print("Downloading RCSB custom report ...")
    url = "https://www.rcsb.org/pdb/rest/customReport.xml?pdbids=*&"
    url += "customReportColumns=%s&" % ",".join(columns)
    url += "format=csv&service=wsfile"
    urllib.request.urlretrieve(url, "pdb-chains.csv")
  print("Reading structures from pdb-chains.csv ...")
  structures = {}
  with open("pdb-chains.csv") as f:
    for row in csv.DictReader(f):
      if any(row[column] == "" for column in columns): continue
      if row["experimentalTechnique"] != "X-RAY DIFFRACTION": continue
      if row["entityMacromoleculeType"] != "Polypeptide(L)": continue
      structureId = row["structureId"]
      if structureId not in structures:
        structures[structureId] = Structure(row)
      structures[structureId].chains.append(Chain(row))
  return list(structures.values())

def get_structure(structure_id):
  for structure in get_structures():
    if structure.id == structure_id:
      return structure
