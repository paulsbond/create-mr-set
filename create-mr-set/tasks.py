import os
import re
import subprocess
import xml.etree.ElementTree as ET


def _run(executable, args=[], stdin=[], stdout=None, stderr=None):
  pstdin = subprocess.PIPE if len(stdin) > 0 else None
  pstdout = None if stdout is None else open(stdout, "w")
  pstderr = None if stderr is None else open(stderr, "w")
  command = [executable] + args
  p = subprocess.Popen(command,
    stdin=pstdin, stdout=pstdout, stderr=pstderr, encoding="utf8")
  if pstdin == subprocess.PIPE:
    for line in stdin:
      p.stdin.write(line + "\n")
    p.stdin.close()
  p.wait()


def mr(hklin, xyzin, identity, prefix, copies, atom_counts):
  """Perform molecular replacement with PHASER"""
  xyzout = "%s.1.pdb" % prefix
  stdout = "%s.log" % prefix
  stderr = "%s.err" % prefix
  keywords = [
    "MODE MR_AUTO",
    "HKLIN %s" % hklin,
    "ENSEMBLE model PDBFILE %s IDENTITY %s" % (xyzin, identity),
    "SEARCH ENSEMBLE model NUM %d" % copies,
    "ROOT %s" % prefix,
    "PURGE ROT NUMBER 1",
    "PURGE TRA NUMBER 1",
    "PURGE RNP NUMBER 1",
    "JOBS 1",
  ]
  for atom in atom_counts:
    keywords.append("COMPOSITION ATOM %-2s NUMBER %d" % (atom, atom_counts[atom]))
  _run("phaser", stdin=keywords, stdout=stdout, stderr=stderr)
  if not os.path.exists(xyzout):
    with open(stdout) as f: log = f.read()
    if not "EXIT STATUS:" in log:
      return { "error": "Early termination" }
    elif "EXIT STATUS: SUCCESS" in log:
      return { "error": "No solution found" }
    elif "INPUT ERROR: No scattering in coordinate file" in log:
      return { "error": "No scattering in input coordinates" }
    elif "INPUT ERROR: Structure Factors of Models" in log:
      return { "error": "Bad ensemble given as input" }
    else:
      return { "error": "No coordinates produced" }
  with open(xyzout) as pdb_file:
    for line in pdb_file:
      if line[:26] == "REMARK Log-Likelihood Gain":
        return { "llg": float(line.split()[-1]) }


def refine(hklin, xyzin, prefix, cycles=10):
  """Refine a structure with REFMAC"""
  hklout = "%s.mtz" % prefix
  xyzout = "%s.pdb" % prefix
  xmlout = "%s.xml" % prefix
  stdout = "%s.log" % prefix
  stderr = "%s.err" % prefix
  _run("refmac5", [
    "HKLIN", hklin,
    "XYZIN", xyzin,
    "HKLOUT", hklout,
    "XYZOUT", xyzout,
    "XMLOUT", xmlout,
  ], [
    "NCYCLES %d" % cycles,
    "PHOUT",
    "END"
  ], stdout=stdout, stderr=stderr)
  with open(stderr) as f:
    for line in f:
      if "Refmac:  " in line:
        return { "error": line }
  for path in (hklout, xyzout, xmlout):
    if not os.path.exists(path):
      return { "error": "Output file missing: %s" % path }
  root = ET.parse(xmlout).getroot()
  rworks = list(root.iter("r_factor"))
  rfrees = list(root.iter("r_free"))
  return {
    "final_rfree": float(rfrees[-1].text),
    "final_rwork": float(rworks[-1].text),
    "initial_rfree": float(rfrees[0].text),
    "initial_rwork": float(rworks[0].text),
  }


def superpose(xyzin1, chain1, xyzin2, chain2, prefix):
  """Superpose one chain over another with GESAMT"""
  xyzout = "%s.pdb" % prefix
  seqout = "%s.seq" % prefix
  stdout = "%s.log" % prefix
  stderr = "%s.err" % prefix
  _run("gesamt", [
    xyzin1, "-s", "//%s" % chain1,
    xyzin2, "-s", "//%s" % chain2,
    "-o", xyzout,
    "-a", seqout,
  ], stdout=stdout, stderr=stderr)
  with open(stdout) as f: log = f.read()
  if "DISSIMILAR and cannot be reasonably aligned" in log:
    return { "error": "GESAMT: Query and target are too dissimilar" }
  for path in (xyzout, seqout):
    if not os.path.exists(path):
      return { "error": "Output file missing: %s" % path }
  result = {}
  match = re.search(r" Q-score          : (\d+\.\d+)", log)
  result["qscore"] = float(match.group(1))
  match = re.search(r" RMSD             : (\d+\.\d+)", log)
  result["rmsd"] = float(match.group(1))
  match = re.search(r" Aligned residues : (\d+)", log)
  result["length"] = int(match.group(1))
  match = re.search(r" Sequence Id:     : (\d+\.\d+)", log)
  result["seqid"] = float(match.group(1))
  return result
