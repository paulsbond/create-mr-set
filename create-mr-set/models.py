import json
import os
import shutil

class _Base:
  def __init__(self, directory):
    self.global_id = directory.replace("/", "_")
    self.directory = directory
    self.metadata = {}
    self.jobs = {}
    os.makedirs(self.directory)

  def path(self, filename):
    return os.path.join(self.directory, filename)

  def add_metadata(self, key, value):
    self.metadata[key] = value
    path = os.path.join(self.directory, "metadata.json")
    with open(path, "w") as f:
      json.dump(self.metadata, f, sort_keys=True, indent=2)

  def remove_directory(self):
    shutil.rmtree(self.directory)

class Structure(_Base):
  def __init__(self, structure_id):
    super().__init__(os.path.join("structures", structure_id))
    self.id = structure_id
    self.chains = {}

class Chain(_Base):
  def __init__(self, chain_id, structure):
    super().__init__(os.path.join(structure.directory, "chains", chain_id))
    self.id = chain_id
    self.homologues = {}
    self.structure = structure
    structure.chains[self.id] = self

class Homologue(_Base):
  def __init__(self, homologue_id, chain):
    super().__init__(os.path.join(chain.directory, "homologues", homologue_id))
    self.id = homologue_id
    self.chain = chain
    chain.homologues[self.id] = self
