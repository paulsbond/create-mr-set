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

def parallel(title, func, structures, processes=None):
  progress_bar = ProgressBar(title, len(structures))
  def callback(s):
    structures[s.id] = s
    progress_bar.increment()
  pool = multiprocessing.Pool(processes)
  for s in structures.values():
    pool.apply_async(func, args=(s,), callback=callback)
  pool.close()
  pool.join()
  progress_bar.finish()

def remove_errors(structures):
  for structureId, structure in list(structures.items()):
    for jobId in structure.jobs:
      if "error" in structure.jobs[jobId]:
        message = "%s: %s" % (jobId, structure.jobs[jobId]["error"])
        structure.add_metadata("error", message)
        del structures[structureId]
  if len(structures) < 1:
    sys.exit("No structures left after removing errors!")
