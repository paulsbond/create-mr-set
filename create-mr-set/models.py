import json
import os

class _Base:
  def __init__(self, directory):
    self.directory = directory
    self.metadata = {}
    self.jobs = {}

  def make_directory(self):
    os.makedirs(self.directory, exist_ok=True)

  def path(self, filename):
    return os.path.join(self.directory, filename)

  def add_metadata(self, key, value):
    self.metadata[key] = value
    path = os.path.join(self.directory, "metadata.json")
    with open(path, "w") as f:
      json.dump(self.metadata, f, sort_keys=True, indent=2)

class Structure(_Base):
  def __init__(self, structure_id):
    super(os.path.join("structures", structure_id))
    self.id = structure_id
    self.chains = {}

class Chain(_Base):
  def __init__(self, chain_id, structure):
    super(os.path.join(structure.directory, "chains", chain_id))
    self.id = chain_id
    self.structure = structure
    self.homologues = {}

class Homologue(_Base):
  def __init__(self, homologue_id, chain):
    super(os.path.join(chain.directory, "homologues", homologue_id))
    self.id = homologue_id
