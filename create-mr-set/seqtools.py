import Bio.Seq
import Bio.SeqIO

def uppercase(seqin, seqout):
  """Transform sequences to uppercase"""
  records = list(Bio.SeqIO.parse(seqin, "fasta"))
  for record in records:
    record.seq = Bio.Seq.Seq(str(record.seq).upper())
  Bio.SeqIO.write(records, seqout, "fasta")
