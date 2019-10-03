import sys

class ProgressBar:
  def __init__(self, task, total, milestones=5):
    self.done = 0
    self.total = total
    step = int(total / milestones)
    self.milestones = {0 + i * step for i in range(milestones)}
    self.active = False
    print("%s..." % task)
    if sys.stdout.isatty():
      self.draw()

  def increment(self):
    self.done += 1
    if sys.stdout.isatty() or self.done in self.milestones:
      self.draw()

  def draw(self):
    hashes = "#" * round(self.done / self.total * 60)
    dashes = "-" * (60 - len(hashes))
    bar = "Progress |%s%s| %d/%d" % (hashes, dashes, self.done, self.total)
    if sys.stdout.isatty():
      if self.active:
        bar = "\r" + bar
      print(bar, end="")
      self.active = True
    else:
      print(bar)

  def finish(self):
    self.draw()
    print("\n" if sys.stdout.isatty() else "")
    self.active = False
