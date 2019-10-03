import subprocess
import os
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

def refine(hklin, xyzin, prefix, cycles=10):
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
