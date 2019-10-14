import csv
import json
import os
import urllib.request

class _Structure:
  def __init__(self, row):
    self.id = row["structureId"]
    self.resolution = float(row["resolution"])
    self.chains = {}
    self.jobs = {}
    self.metadata = {
      "reported_resolution": self.resolution,
      "reported_rwork": float(row["rWork"]),
      "reported_rfree": float(row["rFree"]),
    }

  def add_metadata(self, key, value):
    self.metadata[key] = value
    path = os.path.join(self.id, "metadata.json")
    with open(path, "w") as f:
      json.dump(self.metadata, f, sort_keys=True, indent=2)

class _Chain:
  def __init__(self, row):
    self.id = row["chainId"]
    self.cluster95 = row["clusterNumber95"]
    self.cluster90 = row["clusterNumber90"]
    self.cluster70 = row["clusterNumber70"]
    self.cluster50 = row["clusterNumber50"]
    self.cluster40 = row["clusterNumber40"]
    self.cluster30 = row["clusterNumber30"]
    self.jobs = {}

def download_custom_report(columns, path):
  print("Downloading custom report ...")
  url = "https://www.rcsb.org/pdb/rest/customReport.xml?pdbids=*&"
  url += "customReportColumns=%s&" % ",".join(columns)
  url += "format=csv&service=wsfile"
  urllib.request.urlretrieve(url, path)

def structures():
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
  if not os.path.exists("pdb-chains.csv"):
    download_custom_report(columns, "pdb-chains.csv")
  print("Reading structures from pdb-chains.csv ...")
  structures = {}
  with open("pdb-chains.csv") as f:
    for row in csv.DictReader(f):
      if any(row[column] == "" for column in columns): continue
      if row["experimentalTechnique"] != "X-RAY DIFFRACTION": continue
      if row["entityMacromoleculeType"] != "Polypeptide(L)": continue
      structureId = row["structureId"]
      chainId = row["chainId"]
      if structureId not in structures:
        structures[structureId] = _Structure(row)
      structures[structureId].chains[chainId] = _Chain(row)
  return structures
