import Bio.Seq
import Bio.SeqIO
import multiprocessing
import subprocess
import sys

class ProgressBar:
  def __init__(self, task, total, milestones=5):
    self.done = 0
    self.total = total
    step = int(total / milestones)
    self.milestones = {0 + i * step for i in range(milestones)}
    self.active = False
    print("%s ..." % task)
    if sys.stdout.isatty():
      self.draw()

  def increment(self):
    self.done += 1
    if sys.stdout.isatty() or self.done in self.milestones:
      self.draw()

  def draw(self):
    hashes = "#" * round(self.done / self.total * 60)
    dashes = "-" * (60 - len(hashes))
    bar = "|%s%s| %d/%d" % (hashes, dashes, self.done, self.total)
    if sys.stdout.isatty():
      if self.active:
        bar = "\r" + bar
      print(bar, end="")
      self.active = True
    else:
      print(bar)

  def finish(self):
    self.draw()
    if sys.stdout.isatty():
      print("")
    self.active = False

def run(executable, args=[], stdin=[], stdout=None, stderr=None):
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

def parallel(title, func, dictionary, processes=None):
  progress_bar = ProgressBar(title, len(dictionary))
  def callback(item):
    key, value = item
    dictionary[key] = value
    progress_bar.increment()
  pool = multiprocessing.Pool(processes)
  for key, value in dictionary.items():
    pool.apply_async(func, args=(key, value), callback=callback)
  pool.close()
  pool.join()
  progress_bar.finish()

def remove_errors(structures):
  error_counts = {}
  for structureId, structure in list(structures.items()):
    for jobId in structure.jobs:
      if "error" in structure.jobs[jobId]:
        message = "%s: %s" % (jobId, structure.jobs[jobId]["error"])
        if message not in error_counts: error_counts[message] = 0
        error_counts[message] += 1
        structure.add_metadata("error", message)
        del structures[structureId]
  if len(error_counts) > 0:
    print("Removed some structures due to errors:")
    for error in error_counts:
      print("%s (%d removed)" % (error, error_counts[error]))
  if len(structures) < 1:
    sys.exit("No structures left after removing errors!")

def is_semet(pdbin):
  """Check if a PDB format file is a selenomethione derivative"""
  with open(pdbin) as f:
    for line in f:
      if line[:6] == "ATOM  " or line[:6] == "HETATM":
        if line[17:20] == "MET": return False
        if line[17:20] == "MSE": return True
  return False

def uppercase(seqin, seqout):
  """Transform sequences to uppercase"""
  records = list(Bio.SeqIO.parse(seqin, "fasta"))
  for record in records:
    record.seq = Bio.Seq.Seq(str(record.seq).upper())
  Bio.SeqIO.write(records, seqout, "fasta")
