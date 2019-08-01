import csv
import os
import sys
import utils

class Structure:
  def __init__(self, row):
    self.id = row["structureId"].lower()
    self.rWork = float(row["rWork"])
    self.rFree = float(row["rFree"])
    self.resolution = round(float(row["resolution"]), 1)
    self.chains = []
    self.represetativeChainsIdentified = False

    self.cif2mtz_log = "%s/cif2mtz.log" % self.id
    self.cif2mtz_mtz = "%s/cif2mtz.mtz" % self.id
    self.deposited_cif = "%s/deposited.cif" % self.id
    self.deposited_pdb = "%s/deposited.pdb" % self.id
    self.modeltoseq_log = "%s/modeltoseq.log" % self.id
    self.modeltoseq_seq = "%s/modeltoseq.seq" % self.id
    self.protein_pdb = "%s/protein.pdb" % self.id
    self.protein_seq = "%s/protein.seq" % self.id
    self.refmac_log = "%s/refmac.log" % self.id
    self.refmac_mtz = "%s/refmac.mtz" % self.id
    self.refmac_pdb = "%s/refmac.pdb" % self.id
    self.refmac_xml = "%s/refmac.xml" % self.id
    self.validation_report = "%s/validation_report.xml.gz" % self.id

  def make_dir(self):
    if not os.path.exists(self.id): os.mkdir(self.id)

  def download_validation_report(self):
    url = "http://ftp.ebi.ac.uk/pub/databases/pdb/validation_reports/%s/%s/%s_validation.xml.gz" % (self.id[1:3], self.id, self.id)
    return utils.download(url, self.validation_report)

  def download_deposited_pdb(self):
    url = "https://files.rcsb.org/download/%s.pdb" % self.id
    return utils.download(url, self.deposited_pdb)

  def download_deposited_cif(self):
    url = "https://files.rcsb.org/download/%s-sf.cif" % self.id
    return utils.download(url, self.deposited_cif)

  def identify_representative_chains(self):
    representatives = {}
    for chain in self.chains:
      clusterNumber = chain.clusterNumber95
      if clusterNumber in representatives:
        representatives[clusterNumber].represents.append(chain)
      else:
        chain.representative = True
        representatives[clusterNumber] = chain
    self.represetativeChainsIdentified = True

class Chain:
  def __init__(self, structure, row):
    self.structure = structure
    self.id = row["chainId"]
    self.clusterNumber95 = int(row["clusterNumber95"])
    self.clusterNumber40 = int(row["clusterNumber40"])
    self.representative = False
    self.represents = []

    self.mrbump_seq = "%s/search_%s.seq" % (self.structure.id, self.id)
    self.mrbump_log = "%s/search_%s.log" % (self.structure.id, self.id)

def get_structures_from_pdb():
  url = ("https://www.rcsb.org/pdb/rest/customReport.xml?pdbids=*&"
    "customReportColumns=experimentalTechnique,rWork,rFree,resolution,"
    "clusterNumber95,clusterNumber40,entityMacromoleculeType&"
    "format=csv&service=wsfile")
  path = "pdb-chains.csv"
  utils.download(url, path)
  structures = {}
  with open(path) as f:
    reader = csv.DictReader(f)
    for row in reader:
      if (row["experimentalTechnique"] != "X-RAY DIFFRACTION" or
          row["rWork"] == "" or
          row["rFree"] == "" or
          row["resolution"] == "" or
          row["clusterNumber95"] == "" or
          row["clusterNumber40"] == "" or
          row["entityMacromoleculeType"] != "Polypeptide(L)"):
        continue
      structureId = row["structureId"]
      if structureId not in structures:
        structures[structureId] = Structure(row)
      structures[structureId].chains.append(Chain(structures[structureId], row))
  return list(structures.values())
