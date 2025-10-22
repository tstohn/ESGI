import pathlib

# ---- tweakables ----
FA_FILE = "toy_ref.fa"                                  # your FASTA
HDR_FASTQ = "SIGNALseq_with_toyRNA_2.fastq"         # source of headers to reuse
OUT_FASTQ = "toy_reads_74_renamed.fq"              # output
READ_LEN = 74                                      # simulated read length
STEP = 74                                          # slide step (set <READ_LEN for overlap)
QUAL_CHAR = "I"                                     # Phred+33 quality char (I ~ Q40)
# --------------------

def load_fastq_headers(path):
    hdrs = []
    with open(path, "r", encoding="utf-8") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 0:                   # 1,5,9,... (0-indexed: 0,4,8,...)
                hdrs.append(line.rstrip("\r\n"))
    return hdrs

def slide_emit(handle, seq, hdrs, outfh):
    seq = seq.replace("\n","").replace("\r","")
    n = len(seq)
    idx = 0
    pos = 0
    while pos <= max(0, n - 1):
        read = seq[pos:pos+READ_LEN]
        if not read:
            break
        q = QUAL_CHAR * len(read)
        # header: reuse if available; else synthetic
        if idx < len(hdrs):
            header_line = hdrs[idx]                      # keep exact header line
        else:
            header_line = f"@{handle}_{pos+1}-{pos+len(read)}"
        outfh.write(f"{header_line}\n{read}\n+\n{q}\n")
        idx += 1
        pos += STEP
        if pos > n - READ_LEN and pos < n:              # ensure last partial window if STEP skips past tail
            pos = n - READ_LEN
    return idx  # how many headers consumed

# read all headers to reuse
reuse_headers = load_fastq_headers(HDR_FASTQ)

# walk FASTA and emit
fa_text = pathlib.Path(FA_FILE).read_text(encoding="utf-8")
hdr=None; buf=[]; used=0
with open(OUT_FASTQ, "w", encoding="utf-8") as outfh:
    for line in fa_text.splitlines():
        if not line:
            continue
        if line.startswith(">"):
            if hdr is not None:
                used += slide_emit(hdr, "".join(buf), reuse_headers[used:], outfh)
            hdr = line[1:].strip()
            buf = []
        else:
            buf.append(line.strip())
    if hdr is not None:
        used += slide_emit(hdr, "".join(buf), reuse_headers[used:], outfh)

print(f"Done. Wrote {OUT_FASTQ}")