import hashlib
import streamlit as st
import pandas as pd

st.set_page_config(page_title="FASTA QC Reporter", layout="wide")

# ---------------- Constants ----------------
VALID = set("ATGCN")   # allow N (common in industry FASTA)
DNA_ONLY = set("ATGC")

# ---------------- UI: Input ----------------
st.title("🧬 Yash Bakshi’s FASTA QC Reporter")

st.markdown("### Paste FASTA (mobile users)")
fasta_text = st.text_area(
    "Paste FASTA here (recommended for phone users)",
    height=200,
    placeholder=">seq1\nATGCGTATATAGCG...\n>seq2\nTTGCA...",
    key="paste_fasta"
)

st.markdown("### OR upload FASTA file")
uploaded_file = st.file_uploader(
    "Upload FASTA (.fasta / .fa / .txt)",
    type=["fasta", "fa", "txt"],
    key="upload_fasta"
)

if fasta_text.strip():
    fasta_input = fasta_text
elif uploaded_file is not None:
    fasta_input = uploaded_file.getvalue().decode("utf-8", errors="ignore")
else:
    st.info("Paste FASTA text OR upload a file to continue.")
    st.stop()

st.write(
    "QC report: validation, GC%, complexity flags, duplicates, and downloadable CSV."
)

# ---------------- Helpers ----------------
def clean(seq: str) -> str:
    return "".join(seq.upper().split())

def parse_fasta(text: str):
    records = []
    header = None
    seq_parts = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, clean("".join(seq_parts))))
            header = line[1:].strip() or "Unnamed"
            seq_parts = []
        else:
            seq_parts.append(line)
    if header is not None:
        records.append((header, clean("".join(seq_parts))))
    return records

def gc_percent(seq: str) -> float:
    core = [b for b in seq if b in DNA_ONLY]  # ignore N in denominator
    if not core:
        return 0.0
    s = "".join(core)
    return round((s.count("G") + s.count("C")) / len(s) * 100, 2)

def shannon_entropy(seq: str) -> float:
    from math import log2
    core = [b for b in seq if b in DNA_ONLY]
    if not core:
        return 0.0
    s = "".join(core)
    n = len(s)
    freqs = {b: s.count(b) / n for b in "ATGC"}
    ent = -sum(p * log2(p) for p in freqs.values() if p > 0)
    return round(ent, 3)

def longest_homopolymer(seq: str) -> int:
    best = 0
    cur = 0
    prev = ""
    for ch in seq:
        if ch == prev:
            cur += 1
        else:
            prev = ch
            cur = 1
        best = max(best, cur)
    return best

def md5_hash(seq: str) -> str:
    return hashlib.md5(seq.encode("utf-8")).hexdigest()

# ---------------- Sidebar: Thresholds ----------------
with st.sidebar:
    st.header("QC Thresholds")
    min_len = st.number_input("Min length", min_value=1, value=50, step=10)
    max_len = st.number_input("Max length", min_value=1, value=500000, step=1000)
    max_n_pct = st.slider("Max N% allowed", 0.0, 100.0, 5.0, 0.5)
    min_entropy = st.slider("Min entropy (complexity)", 0.0, 2.0, 1.2, 0.05)
    max_homopolymer = st.number_input("Max homopolymer run", min_value=1, value=10, step=1)

# ---------------- Parse ----------------
records = parse_fasta(fasta_input)

if not records:
    st.error("No FASTA records found. Make sure file has lines starting with '>' and sequences below them.")
    st.stop()

# ---------------- QC + Duplicates ----------------
hash_to_names = {}
rows = []

for name, seq in records:
    length = len(seq)
    invalid = sorted(set(seq) - VALID)
    n_count = seq.count("N")
    n_pct = round((n_count / length * 100), 2) if length else 0.0
    gc = gc_percent(seq)
    ent = shannon_entropy(seq)
    homo = longest_homopolymer(seq)
    h = md5_hash(seq)

    hash_to_names.setdefault(h, []).append(name)

    flags = []
    if invalid:
        flags.append("INVALID_CHARS")
    if length < min_len:
        flags.append("TOO_SHORT")
    if length > max_len:
        flags.append("TOO_LONG")
    if n_pct > max_n_pct:
        flags.append("TOO_MANY_N")
    if ent < min_entropy:
        flags.append("LOW_COMPLEXITY")
    if homo > max_homopolymer:
        flags.append("HOMOPOLYMER_RUN")

    status = "PASS" if len(flags) == 0 else "FAIL"

    rows.append({
        "sequence_name": name,
        "length": length,
        "gc_percent": gc,
        "n_count": n_count,
        "n_percent": n_pct,
        "entropy": ent,
        "max_homopolymer": homo,
        "invalid_chars": "".join(invalid),
        "md5": h,
        "status": status,
        "flags": ",".join(flags)
    })

df = pd.DataFrame(rows)

# Mark duplicates
dup_map = {}
for h, names in hash_to_names.items():
    if len(names) > 1:
        for n in names:
            dup_map[n] = "|".join(names)

df["duplicates"] = df["sequence_name"].map(dup_map).fillna("")

# ---------------- Output ----------------
st.subheader("📌 QC Summary")
pass_count = int((df["status"] == "PASS").sum())
fail_count = int((df["status"] == "FAIL").sum())
st.write(f"✅ PASS: **{pass_count}**   |   ❌ FAIL: **{fail_count}**   |   Total: **{len(df)}**")

st.subheader("📄 QC Report Table")
st.dataframe(df.sort_values(["status", "sequence_name"]), use_container_width=True)

st.subheader("📈 Quick Charts")
c1, c2 = st.columns(2)
with c1:
    st.bar_chart(df.set_index("sequence_name")["gc_percent"])
with c2:
    st.bar_chart(df.set_index("sequence_name")["n_percent"])

csv_bytes = df.to_csv(index=False).encode("utf-8")
st.download_button(
    "⬇️ Download qc_report.csv",
    data=csv_bytes,
    file_name="qc_report.csv",
    mime="text/csv"
)

st.caption("Tip: This QC report is designed like a preprocessing step used in real bioinformatics pipelines.")

