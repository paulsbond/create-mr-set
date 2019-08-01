import distutils.spawn
import os
import subprocess
import sys
import xml.etree.ElementTree as ET

# TODO: Remove PIPE from stdout and stderr
def run(command, stdin = [], stdout = None, stderr = None):
  pstdin = subprocess.PIPE if len(stdin) > 0 else None
  pstdout = subprocess.PIPE if stdout is None else open(stdout, 'w')
  pstderr = subprocess.PIPE if stderr is None else open(stderr, 'w')
  command = command.split(' ')
  p = subprocess.Popen(command,
    stdin=pstdin, stdout=pstdout, stderr=pstderr, encoding='utf8')
  if pstdin == subprocess.PIPE:
    for line in stdin:
      p.stdin.write(line + "\n")
    p.stdin.close()
  p.wait()

def ensure_executables_exist():
  executables = [
    "cextractprotein", "cif2mtz", "cmodeltoseq", "gesamt", "mrbump", "pdbcur",
    "pdb_merge", "refmac5", "phaser.sculptor"]
  for executable in executables:
    filename = distutils.spawn.find_executable(executable)
    if filename == None:
      sys.exit("Executable '%s' not found." % executable)

def count_atoms(pdbin):
  command = ["grep", "-c", "^ATOM  ", pdbin]
  p = subprocess.Popen(command, stdout=subprocess.PIPE)
  return int(p.stdout.read())

def cextractprotein(xyzin, xyzout, logfile = None):
  command = "cextractprotein"
  command += " -xyzin " + xyzin
  command += " -xyzout " + xyzout
  command += " -min-length 6"
  command += " -truncate-unk"
  run(command, [], logfile)

def cif2mtz(hklin, hklout, logfile = None):
  command = "cif2mtz hklin %s hklout %s" % (hklin, hklout)
  stdin = ["END"]
  run(command, stdin, logfile)

def cmodeltoseq(pdbin, seqout, logfile = None):
  command = "cmodeltoseq"
  command += " -pdbin " + pdbin
  command += " -seqout " + seqout
  run(command, [], logfile)

def extract_chain(seqin, chain, seqout):
  should_extract = False
  f = open(seqout, "w")
  for line in [l for l in open(seqin)]:
    if line[0] == ">":
      should_extract = line[1:].strip() == chain
    if should_extract:
      f.write(line)

def gesamt(xyzref, refchain, xyzwrk, xyzout, logfile=None):
  command = "gesamt %s -s /*/%s %s -high -o %s" % (xyzref, refchain, xyzwrk, xyzout)
  run(command, [], logfile)

def mrbump(chain_id, structure_id, seqin, logfile):
  command = "mrbump"
  command += " SEQIN " + seqin
  stdin = [
    "JOBID %s" % chain_id,
    "ROOTDIR %s" % structure_id,
    "RLEVEL 95",
    "MRNUM 100",
    "MDLMOLREP False",
    "MDLCHAINSAW False",
    "LITE True",
    "USEENSEM False",
    "END"
  ]
  run(command, stdin, logfile)

def pdbcur(xyzin, xyzout, stdin):
  command = "pdbcur"
  command += " xyzin " + xyzin
  command += " xyzout " + xyzout
  run(command, stdin)

def pdb_merge(xyzin1, xyzin2, xyzout):
  command = "pdb_merge xyzin1 %s xyzin2 %s xyzout %s" % (xyzin1, xyzin2, xyzout)
  stdin = ["NOMERGE", "END"]
  run(command, stdin)

def refmac5(hklin, xyzin, hklout, xyzout, xmlout, cycles, logfile = None):
  command = "refmac5"
  command += " HKLIN " + hklin + " HKLOUT " + hklout
  command += " XYZIN " + xyzin + " XYZOUT " + xyzout + " XMLOUT " + xmlout
  stdin = [
    "NCYCLES %d" % cycles,
    "WEIGHT AUTO",
    "MAKE HYDR ALL",
    "REFI BREF ISOT",
    "MAKE NEWLIGAND EXIT",
    "SCALE TYPE SIMPLE",
    "SOLVENT YES",
    "PHOUT",
    "END"
  ]
  run(command, stdin, logfile)

def refmac_terminated_normally(xmlout):
  if not os.path.exists(xmlout): return False
  root = ET.parse(xmlout).getroot()
  err_level = root.find("err_level").text.strip()
  err_message = root.find("err_message").text.strip()
  return err_level == "0" and err_message == "normal termination"

def sculptor(xyzin, alignin, folder, domain, logfile=None):
  command = "phaser.sculptor --stdin"
  stdin = [
    "input {",
    "model {"
    "file_name = %s" % xyzin,
    "remove_alternate_conformations = True"
    "}",
    "alignment {",
    "file_name = %s" % alignin,
    "target_index = 1",
    "}",
    "}",
    "output {",
    "folder = %s" % folder,
    "root = 'sculptor%s'" % domain["id"],
    "}",
    "macromolecule {",
    "pruning {",
    "use = schwarzenbacher",
    "}",
    "bfactor {",
    "use = asa",
    "}",
    "}"
  ]
  run(command, stdin, logfile)
