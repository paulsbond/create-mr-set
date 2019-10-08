import Bio.Seq
import Bio.SeqIO
import urllib.request
import uuid

def download(pdb_ids, path):
  """Download a FASTA file from RCSB for one or more entries"""
  url = "https://www.rcsb.org/pdb/download/downloadFastaFiles.do"
  url += "?structureIdList=%s" % ",".join(pdb_ids)
  url += "&compressionType=uncompressed"
  urllib.request.urlretrieve(url, path)

def extract_entry(seqin, pdb_id, seqout):
  """Extract the sequences for a single entry from an RCSB FASTA file"""
  pdb_id = pdb_id.upper()
  records = list(Bio.SeqIO.parse(seqin, "fasta"))
  records = [r for r in records if r.id[:4] == pdb_id]
  Bio.SeqIO.write(records, seqout, "fasta")

def uppercase(seqin, seqout):
  """Transform sequences to uppercase"""
  records = list(Bio.SeqIO.parse(seqin, "fasta"))
  for record in records:
    record.seq = Bio.Seq.Seq(str(record.seq).upper())
  Bio.SeqIO.write(records, seqout, "fasta")
