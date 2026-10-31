#!/usr/bin/env python3
"""
read-context.py - which library kit/adapters are still on the reads.

ONT (-p ont, default): looks for the known kit adapters (indels allowed, ~15% edits) in the first 160 nt of a sample of reads, spread
over all files. For barcoded kits it also calls the ONT barcode. The output is a self-contained HTML section for
faster-report (shown in an iframe, like the gel simulation).
Illumina (-p illumina): looks for adapter read-through (TruSeq, Nextera, small RNA) anywhere in the reads and reads
the indexes from the CASAVA 1.8+ read headers.
ONT part ported from `utag context` (https://github.com/angelovangel/utag); stdlib only.

Usage:
    read-context.py reads1.fastq.gz reads2.fastq.gz -o read_context.html [-p ont|illumina] [-m short|long] [-n 5000] [-t 4]
"""

import argparse, gzip, itertools, json, math, os, random, re, sys, time
from collections import Counter
from multiprocessing import Pool

COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s):
    return s.translate(COMP)[::-1]


# --- ONT kit sequences, from barbell src/kits/kits.rs (Dorado barcode_kits), as in utag ---
NB_FRONT = "ATTGCTAAGGTTAA"   # native adapter side, precedes the NB barcode
NB_REAR = "CAGCACCT"         # constant 8 nt directly after the NB barcode
# Ligation Adapter top strand; the final T is the 3' T overhang on the dA tail. Read starts are often miscalled,
# so only the last 22 nt are searched.
LSK_ADAPTER = "CCTGTACTTCGTTCAGTTACGTATTGCT"
LSK_ADAPTER_TAIL = LSK_ADAPTER[-22:]
MAB_TOPO = "CCATATCCGTGTCGCCCTT"              # vaccinia TOPO site after the AB barcode (SQK-MAB114)
MAB_FRONT = "GCTTGGGTGTTTAACC"                 # end of the rapid adapter, directly 5' of the AB / RB barcode
RAPID_FLANK = "GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA"   # transposase flank (SQK-RAD114 / RBK114)
MAB_SEQS = {    # Dorado TP01-TP24 = barbell AB01-AB24
    "AB01": "GCACCTGGAACTTGTGCCTTCCAC", "AB02": "CCGAAATAGGTTATCTGTTGTTGT",
    "AB03": "ATCAATCGCTGGACGATGGATTAG", "AB04": "CCACCCGCTCCTGCCGGTGGGCGT",
    "AB05": "AGACTCTTGGGCTCGCCACGTCCC", "AB06": "TCTGTATCCGGAGACGGGATGGAC",
    "AB07": "TTTCGGATCAATCGACCGCAAACG", "AB08": "ACTCAAACATTCTGTTAGATCGCG",
    "AB09": "AAATGGAACCCGGATATGTTTACT", "AB10": "TAAATCGACCTATGATGAACACAG",
    "AB11": "ACATGTTGGAGTGAAAGTCGGGTA", "AB12": "CCTGGACCACGATCATTGTAACAT",
    "AB13": "TATGGTGGATCTCCCTCTATCTTC", "AB14": "AAGTAAATGGGACGCCCACTCCGA",
    "AB15": "TGTTCGCGGCTTGATCTAATATTA", "AB16": "AGAGAGCTTCCCGGGAGGGTGGTC",
    "AB17": "TTGTGAATATCTGTCACAAACACC", "AB18": "CAATCGTACCAGGGAACATAAAGT",
    "AB19": "CACACCCAAACAATATGGACCCGT", "AB20": "AATAACCACATCCGCCCTCCGCAC",
    "AB21": "TCCTAATAATGTGTAGATCGGTCC", "AB22": "AGTCGATGGAACAAGAGAAGTTAT",
    "AB23": "AAACTCACTGTATGTCGTTTCTAT", "AB24": "TGACATCACTGATCGAGGAAGATC",
}
NB_SEQS = {
    "NB01": "CACAAAGACACCGACAACTTTCTT", "NB02": "ACAGACGACTACAAACGGAATCGA",
    "NB03": "CCTGGTAACTGGGACACAAGACTC", "NB04": "TAGGGAAACACGATAGAATCCGAA",
    "NB05": "AAGGTTACACAAACCCTGGACAAG", "NB06": "GACTACTTTCTGCCTTTGCGAGAA",
    "NB07": "AAGGATTCATTCCCACGGTAACAC", "NB08": "ACGTAACTTGGTTTGTTCCCTGAA",
    "NB09": "AACCAAGACTCGCTGTGCCTAGTT", "NB10": "GAGAGGACAAAGGTTTCAACGCTT",
    "NB11": "TCCATTCCCTCCGATAGATGAAAC", "NB12": "TCCGATTCTGCTTCTTTCTACCTG",
    "NB13": "AGAACGACTTCCATACTCGTGTGA", "NB14": "AACGAGTCTCTTGGGACCCATAGA",
    "NB15": "AGGTCTACCTCGCTAACACCACTG", "NB16": "CGTCAACTGACAGTGGTTCGTACT",
    "NB17": "ACCCTCCAGGAAAGTACCTCTGAT", "NB18": "CCAAACCCAACAACCTAGATAGGC",
    "NB19": "GTTCCTCGTGCAGTGTCAAGAGAT", "NB20": "TTGCGTCCTGTTACGAGAACTCAT",
    "NB21": "GAGCCTCTCATTGTCCGTTCTCTA", "NB22": "ACCACTGCCATGTATCAAAGTACG",
    "NB23": "CTTACTACCCAGTGAACCTCCTCG", "NB24": "GCATAGTTCTGCATGATGGGTTAG",
    "NB25": "GTAAGTTGGGTATGCAACGCAATG", "NB26": "CATACAGCGACTACGCATTCTCAT",
    "NB27": "CGACGGTTAGATTCACCTCTTACA", "NB28": "TGAAACCTAAGAAGGCACCGTATC",
    "NB29": "CTAGACACCTTGGGTTGACAGACC", "NB30": "TCAGTGAGGATCTACTTCGACCCA",
    "NB31": "TGCGTACAGCAATCAGTTACATTG", "NB32": "CCAGTAGAAGTCCGACAACGTCAT",
    "NB33": "CAGACTTGGTACGGTTGGGTAACT", "NB34": "GGACGAAGAACTCAAGTCAAAGGC",
    "NB35": "CTACTTACGAAGCTGAGGGACTGC", "NB36": "ATGTCCCAGTTAGAGGAGGAAACA",
    "NB37": "GCTTGCGATTGATGCTTAGTATCA", "NB38": "ACCACAGGAGGACGATACAGAGAA",
    "NB39": "CCACAGTGTCAACTAGAGCCTCTC", "NB40": "TAGTTTGGATGACCAAGGATAGCC",
    "NB41": "GGAGTTCGTCCAGAGAAGTACACG", "NB42": "CTACGTGTAAGGCATACCTGCCAG",
    "NB43": "CTTTCGTTGTTGACTCGACGGTAG", "NB44": "AGTAGAAAGGGTTCCTTCCCACTC",
    "NB45": "GATCCAACAGAGATGCCTTCAGTG", "NB46": "GCTGTGTTCCACTTCATTCTCCTG",
    "NB47": "GTGCAACTTTCCCACAGGTAGTTC", "NB48": "CATCTGGAACGTGGTACACCTGTA",
    "NB49": "ACTGGTGCAGCTTTGAACATCTAG", "NB50": "ATGGACTTTGGTAACTTCCTGCGT",
    "NB51": "GTTGAATGAGCCTACTGGGTCCTC", "NB52": "TGAGAGACAAGATTGTTCGTGGAC",
    "NB53": "AGATTCAGACCGTCTCATGCAAAG", "NB54": "CAAGAGCTTTGACTAAGGAGCATG",
    "NB55": "TGGAAGATGAGACCCTGATCTACG", "NB56": "TCACTACTCAACAGGTGGCATGAA",
    "NB57": "GCTAGGTCAATCTCCTTCGGAAGT", "NB58": "CAGGTTACTCCTCCGTGAGTCTGA",
    "NB59": "TCAATCAAGAAGGGAAAGCAAGGT", "NB60": "CATGTTCAACCAAGGCTTCTATGG",
    "NB61": "AGAGGGTACTATGTGCCTCAGCAC", "NB62": "CACCCACACTTACTTCAGGACGTA",
    "NB63": "TTCTGAAGTTCCTGGGTCTTGAAC", "NB64": "GACAGACACCGTTCATCGACTTTC",
    "NB65": "TTCTCAGTCTTCCTCCAGACAAGG", "NB66": "CCGATCCTTGTGGCTTCTAACTTC",
    "NB67": "GTTTGTCATACTCGTGTGCTCACC", "NB68": "GAATCTAAGCAAACACGAAGGTGG",
    "NB69": "TACAGTCCGAGCCTCATGTGATCT", "NB70": "ACCGAGATCCTACGAATGGAGTGT",
    "NB71": "CCTGGGAGCATCAGGTAGTAACAG", "NB72": "TAGCTGACTGTCTTCCATACCGAC",
    "NB73": "AAGAAACAGGATGACAGAACCCTC", "NB74": "TACAAGCATCCCAACACTTCCACT",
    "NB75": "GACCATTGTGATGAACCCTGTTGT", "NB76": "ATGCTTGTTACATCAACCCTGGAC",
    "NB77": "CGACCTGTTTCTCAGGGATACAAC", "NB78": "AACAACCGAACCTTTGAATCAGAA",
    "NB79": "TCTCGGAGATAGTTCTCACTGCTG", "NB80": "CGGATGAACATAGGATAGCGATTC",
    "NB81": "CCTCATCTTGTGAAGTTGTTTCGG", "NB82": "ACGGTATGTCGAGTTCCAGGACTA",
    "NB83": "TGGCTTGATCTAGGTAAGGTCGAA", "NB84": "GTAGTGGACCTAGAACCTGTGCCA",
    "NB85": "AACGGAGGAGTTAGTTGGATGATC", "NB86": "AGGTGATCCCAACAAGCGTAAGTA",
    "NB87": "TACATGCTCCTGTTGTTAGGGAGG", "NB88": "TCTTCTACTACCGATCCGAAGCAG",
    "NB89": "ACAGCATCAATGTTTGGCTAGTTG", "NB90": "GATGTAGAGGGTACGGTTTGAGGC",
    "NB91": "GGCTCCATAGGAACTCACGCTACT", "NB92": "TTGTGAGTGGAAAGATACAGGACC",
    "NB93": "AGTTTCCATCACTTCAGACTTGGG", "NB94": "GATTGTCCTCAAACTGCCACCTAC",
    "NB95": "CCTGTCTGGAAGAAGAATGGACTT", "NB96": "CTGAACGGTCATAGAGTCCACCAT",
}
# RBK114 (barbell RBK01-RBK96): RB01-RB12 are the reverse complements of NB01-NB12, RB13-RB96 the NB sequences,
# except six with their own sequence
RB_SEQS = {"RB" + k[2:]: rc(v) if int(k[2:]) <= 12 else v for k, v in NB_SEQS.items()}
RB_SEQS.update({"RB26": "ACTATGCCTTTCCGTGAAACAGTT", "RB39": "TCTGCCACACACTCGTAAGTCCTT",
                "RB40": "GTCGATACTGGACCTATCCCTTGG", "RB48": "GAGTCCGTGACAACTTCTGAAAGC",
                "RB54": "GGGTGCCAACTACATACCAAACCT", "RB60": "GAACCCTACTTTGGACAGACACCT"})

KIT_WINDOW = 160                         # read-start nt checked for the kit adapters
MIN_PER_FILE = 50                        # reads checked per file at least (-n is spread over the files)
BARCODE_PER_FILE, BARCODE_MAX = 30, 3000 # kit reads per file checked for the ONT barcode (96 alignments each)
ALIGN_MAX = 2000                         # kit reads in the base-frequency logo
KIT_ADAPTERS = {"la": LSK_ADAPTER_TAIL, "nb_front": NB_FRONT, "nb_rear": NB_REAR,
                "rapid_flank": RAPID_FLANK, "topo": MAB_TOPO}
# kit -> (layout, barcode set or None); checked in this order, the first that fits a read wins
KITS = {"SQK-NBD114": ("ligation adapter + NB front + NB barcode + NB rear", NB_SEQS),
        "SQK-MAB114": ("(rapid adapter + AB barcode +) TOPO site", MAB_SEQS),
        "SQK-RAD114": ("rapid adapter + transposase flank, no barcode", None),
        "SQK-RBK114": ("(rapid adapter + RB barcode +) transposase flank", RB_SEQS),
        "SQK-LSK114": ("ligation adapter only", None)}
# what sits 5' of the insert, 5' -> 3': [label, sequence or None (= the barcode), length]
KIT_LAYOUT = {"SQK-NBD114": [["NB front", NB_FRONT, len(NB_FRONT)], ["NB barcode", None, 24], ["NB rear", NB_REAR, len(NB_REAR)]],
              "SQK-MAB114": [["rapid adapter", MAB_FRONT, len(MAB_FRONT)], ["AB barcode", None, 24], ["TOPO site", MAB_TOPO, len(MAB_TOPO)]],
              "SQK-LSK114": [["Ligation Adapter", LSK_ADAPTER[:-1], len(LSK_ADAPTER) - 1], ["T", "T", 1]],
              "SQK-RBK114": [["rapid adapter", MAB_FRONT, len(MAB_FRONT)], ["RB barcode", None, 24],
                             ["transposase flank", RAPID_FLANK, len(RAPID_FLANK)]],
              "SQK-RAD114": [["rapid adapter", MAB_FRONT, len(MAB_FRONT)], ["transposase flank", RAPID_FLANK, len(RAPID_FLANK)]]}
EXT_RE = re.compile(r"\.(fastq|fq|fasta|fa)(\.gz)?$")


# ---------------------------------------------------------------- alignment --
def suffix_edits(adapter, up):
    """d[k] = fewest edits (indels allowed) to align adapter[-k:] so that it ends exactly where `up` ends; it may
    start anywhere in `up`. One DP over the reversed strings gives every suffix length at once."""
    q, t = adapter[::-1], up[::-1][:len(adapter) + 12]
    row = list(range(len(t) + 1))
    d = [0]
    for x in q:
        new = [row[0] + 1]
        for j, y in enumerate(t, 1):
            new.append(min(row[j] + 1, new[j - 1] + 1, row[j - 1] + (x != y)))
        row = new
        d.append(min(row))
    return d


def seq_start(p, text, k):
    """(start, edits) of the whole sequence p in text within k edits, or None. DP on the reversed strings: the end
    of p is free, so the best row-end gives where it starts."""
    q, t = p[::-1], text[::-1]
    prev = [0] * (len(t) + 1)
    for i, a in enumerate(q, 1):
        cur = [i]
        for j, b in enumerate(t, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (b != a)))
        prev = cur
    e = min(prev)
    return (len(t) - prev.index(e), e) if e <= k else None


def locate(p, text, k):
    """(start, end, edits) of p in text within k edits (indels allowed), or None. Exact seeds (k+1 pieces: one is
    exact, pigeonhole) pick the candidates; a semi-global DP around each gives the end, seq_start the start."""
    b = [i * len(p) // (k + 1) for i in range(k + 2)]
    cands, best = sorted({m.start() - b[i] for i in range(k + 1)
                          for m in re.finditer("(?=" + p[b[i]:b[i + 1]] + ")", text)}), None
    done = -1
    for c in cands:
        if c <= done:
            continue
        lo, hi = max(0, c - k), min(len(text), c + len(p) + k)
        done = c + k
        t = text[lo:hi]
        row = [0] * (len(t) + 1)
        for i, x in enumerate(p, 1):
            new = [i]
            for j, y in enumerate(t, 1):
                new.append(min(row[j] + 1, new[j - 1] + 1, row[j - 1] + (x != y)))
            row = new
        e = min(row)
        if e <= k and (best is None or e < best[2]):
            end = lo + row.index(e)
            best = (lo + seq_start(p, t[:end - lo], k)[0], end, e)
    return best


def present_edits(n):
    """edits allowed for a known adapter to count as present (~15%, at least 1)"""
    return max(1, round(0.15 * n))


def read_kit(x):
    """(kit, read up to the end of its ONT barcode, read up to the insert) for the start of read x;
    (None, None, None) if no kit fits"""
    f = {key: locate(seq, x, present_edits(len(seq))) for key, seq in KIT_ADAPTERS.items()}
    if f["nb_front"] and f["nb_rear"] and 18 <= f["nb_rear"][0] - f["nb_front"][1] <= 30:   # a 24 nt barcode fits
        return "SQK-NBD114", x[:f["nb_rear"][0]], x[:f["nb_rear"][1]]
    if f["topo"]:                               # distinctive; the rapid adapter before the AB barcode is often garbled
        return "SQK-MAB114", x[:f["topo"][0]], x[:f["topo"][1]]
    if f["rapid_flank"]:                        # rapid adapter end right before the flank (RAD) or ~24 nt before (RBK)
        r = x[:f["rapid_flank"][0]]
        rad = suffix_edits(MAB_FRONT, r)[-1]
        rbk = min(suffix_edits(MAB_FRONT, r[:-n])[-1] for n in range(20, 29))
        return ("SQK-RAD114", None, x[:f["rapid_flank"][1]]) if rad <= rbk else ("SQK-RBK114", r, x[:f["rapid_flank"][1]])
    if f["la"]:
        return "SQK-LSK114", None, x[:f["la"][1]]
    return None, None, None


def best_barcode(seqs, regions):
    """Which ONT barcode ends where each region ends: Counter of barcode -> reads ('?' = unclear)"""
    c = Counter()
    for region in regions:
        scores = sorted((suffix_edits(s, region)[-1], b) for b, s in seqs.items())
        best, second = scores[0], scores[1]
        c[best[1] if best[0] <= 6 and second[0] - best[0] >= 2 else "?"] += 1
    return c


def base_profile(L):
    """per position of the equal-length strings L ('-' = no base): ([A, C, G, T] share of all strings, consensus
    base or '.', % of all strings with the consensus base)"""
    cons, supp, logo = [], [], []
    for pos in range(len(L[0]) if L else 0):
        c = Counter(x[pos] for x in L if x[pos] != "-")
        logo.append([round(c[b] / len(L), 4) for b in "ACGT"])
        n = sum(c.values())
        b, k = c.most_common(1)[0] if c else ("-", 0)
        cons.append(b if n and k / len(L) >= 0.5 else ".")
        supp.append(round(100 * k / len(L), 1) if n else 0)
    return logo, "".join(cons), supp


# ---------------------------------------------------------------- reading --
def read_starts(path):
    """first KIT_WINDOW nt of each read in a fastq/fasta(.gz) file"""
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        head = fh.readline()
        if head.startswith(">"):                # fasta: a record may span several lines
            seq = ""
            for ln in fh:
                if ln.startswith(">"):
                    yield seq[:KIT_WINDOW].upper()
                    seq = ""
                elif len(seq) < KIT_WINDOW:
                    seq += ln.strip()
            yield seq[:KIT_WINDOW].upper()
            return
        while head:
            yield fh.readline().strip()[:KIT_WINDOW].upper()
            fh.readline(); fh.readline()
            head = fh.readline()


def scan_file(job):
    """kit calls for the first `take` reads of a file: [(kit, barcode region, insert region)]"""
    path, take = job
    return [read_kit(x) for x in itertools.islice(read_starts(path), take) if x]


def barcode_file(job):
    seqs, regions = job
    return best_barcode(seqs, regions)


def label(path):
    return EXT_RE.sub("", os.path.basename(path))


# ---------------------------------------------------------------- main --
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("reads", nargs="+", help="fastq/fasta files (.gz ok)")
    ap.add_argument("-o", "--out", default="read_context.html", help="HTML output (default read_context.html)")
    ap.add_argument("-n", "--n-reads", type=int, default=5000,
                    help=f"reads checked in total, spread over the files (at least {MIN_PER_FILE} per file; default 5000)")
    ap.add_argument("-t", "--threads", type=int, default=4, help="processes (default 4)")
    ap.add_argument("-m", "--mode", choices=["short", "long"], default="long",
                    help="short: only the kit (and ONT barcodes) found; long: also the tables, logo and per-file view (default)")
    ap.add_argument("-p", "--platform", choices=["ont", "illumina"], default="ont",
                    help="ont: kit adapters/barcodes at the read starts (default); illumina: adapter read-through and "
                         "the indexes in the read headers")
    ap.add_argument("--json", help="also write the data as JSON")
    a = ap.parse_args()

    files = sorted(a.reads, key=lambda p: [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", label(p))])
    per = max(MIN_PER_FILE, math.ceil(a.n_reads / len(files)))
    if a.platform == "illumina":
        return illumina(a, files, per)
    with Pool(max(1, min(a.threads, len(files)))) as pool:
        calls = pool.map(scan_file, [(f, per) for f in files])

        c = Counter(k for cl in calls for k, _, _ in cl)
        n = max(1, sum(c.values()))
        R = {"mode": a.mode, "started": time.strftime("%Y-%m-%d %H:%M"), "n": sum(c.values()), "window": KIT_WINDOW, "files": len(files),
             "per_file": per, "kits": [[k, lay, c[k] / n] for k, (lay, _) in KITS.items()] + [["none", "", c[None] / n]]}
        kit, m = max(((k, c[k]) for k in KITS), key=lambda kv: kv[1])
        if m >= 0.2 * n:
            R["kit"] = kit
            seqs = KITS[kit][1]
            bcs = [None] * len(files)
            if seqs:                            # which ONT barcodes: a few kit reads per file
                cap = max(5, min(BARCODE_PER_FILE, BARCODE_MAX // len(files)))
                bcs = pool.map(barcode_file, [(seqs, [r for k, r, _ in cl if k == kit][:cap]) for cl in calls])
                bc = sum(bcs, Counter())
                nc = max(1, sum(bc.values()))
                R.update(barcodes=[[b, v / nc, seqs[b]] for b, v in bc.most_common() if b != "?"],
                         checked=sum(bc.values()), unclear=bc["?"] / nc, set=list(seqs.items()))
            ups = [u for cl in calls for k, _, u in cl if k == kit]
            random.Random(1).shuffle(ups)       # not only the first files
            lay = KIT_LAYOUT[kit]
            U = min(KIT_WINDOW, sum(x for _, _, x in lay) + 10)
            logo, cons, supp = base_profile([u[-U:].rjust(U, "-") for u in ups[:ALIGN_MAX]])
            R["align"] = {"logo": logo, "consensus": cons, "support": supp, "layout": lay, "reads": min(len(ups), ALIGN_MAX)}

    R["file_rows"] = []
    for f, cl, i in zip(files, calls, range(len(files))):
        fc = Counter(k for k, _, _ in cl)
        fn = max(1, len(cl))
        row = {"file": label(f), "n": len(cl), "shares": {(k or "none"): v / fn for k, v in fc.items()}}
        if R.get("barcodes") is not None and bcs[i]:
            b = bcs[i]
            bn = sum(b.values())
            row["bc"] = [[k, v / bn] for k, v in b.most_common() if k != "?"][:3]
            row["bc_n"] = bn
        R["file_rows"].append(row)

    if "kit" not in R:
        R["verdict"] = ("no known ONT adapter at the read starts: the reads are trimmed")
    elif R["kit"] == "SQK-RBK114":
        R["verdict"] = "rapid kit: the transposase adapter and barcodes are still on the reads (untrimmed); inserts start at random positions"
    elif R["kit"] == "SQK-RAD114":              # no barcode in this kit
        R["verdict"] = "rapid kit: the transposase adapter is still on the reads (untrimmed); inserts start at random positions"
    elif R["kit"] in ("SQK-NBD114", "SQK-MAB114"):
        R["verdict"] = "the kit adapter and barcodes are still on the reads (untrimmed)"
    else:                                       # SQK-LSK114: no barcode
        R["verdict"] = "the kit adapter is still on the reads (untrimmed)"

    print(f"read context: {R['n']} reads from {len(files)} files; " +
          ", ".join(f"{k} {100 * f:.1f}%" for k, _, f in R["kits"]) + f" -> {R.get('kit', 'none')}", file=sys.stderr)
    write(a, R, ONT_BODY, ONT_JS)


def write(a, R, body, js):
    if a.json:
        with open(a.json, "w") as fh:
            json.dump(R, fh)
    with open(a.out, "w") as fh:
        fh.write(page(body, js).replace("__DATA__", json.dumps(R, separators=(",", ":")).replace("</", "<\\/")))


# ---------------------------------------------------------------- illumina --
# Illumina reads start with the insert; the adapter shows up at the 3' end when the insert is shorter than the read
# (read-through). family -> (library, adapter start shared by read 1 and read 2, edits allowed)
ILLUMINA_ADAPTERS = {"TruSeq": ("TruSeq / ligation (Illumina DNA PCR-Free, Stranded mRNA, TruSeq)", "AGATCGGAAGAGC", 1),
                     "Nextera": ("Nextera / tagmentation (Illumina DNA Prep, Nextera XT)", "CTGTCTCTTATACACATCT", 2),
                     "small RNA": ("TruSeq Small RNA 3' adapter", "TGGAATTCTCGG", 1)}
ILLUMINA_FOUND = 0.005                   # share of reads with an adapter for it to count as found
INDEX_RE = re.compile(r"\s\d:[YN]:\d+:([ACGTN]+(?:\+[ACGTN]+)?)\s*$")   # CASAVA 1.8+: 1:N:0:ATCACG+TTAGGC


def read_records(path):
    """(header, sequence) of each read in a fastq/fasta(.gz) file"""
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        head = fh.readline()
        if head.startswith(">"):
            seq = []
            for ln in fh:
                if ln.startswith(">"):
                    yield head, "".join(seq).upper()
                    head, seq = ln, []
                else:
                    seq.append(ln.strip())
            yield head, "".join(seq).upper()
            return
        while head:
            yield head, fh.readline().strip().upper()
            fh.readline(); fh.readline()
            head = fh.readline()


def read_through(x):
    """(adapter family, read position where it starts) of the first known adapter in read x, or (None, None)"""
    hits = [(h[0], fam) for fam, (_, sq, k) in ILLUMINA_ADAPTERS.items() for h in [locate(sq, x, k)] if h]
    return (min(hits)[1], min(hits)[0]) if hits else (None, None)


def scan_illumina(job):
    """adapter calls [(family, position, read length)] and index Counter for the first `take` reads of a file"""
    path, take = job
    calls, idx = [], Counter()
    for head, seq in itertools.islice(read_records(path), take):
        if not seq:
            continue
        m = INDEX_RE.search(head)
        idx[m.group(1) if m else None] += 1
        calls.append((*read_through(seq), len(seq)))
    return calls, idx


def illumina(a, files, per):
    with Pool(max(1, min(a.threads, len(files)))) as pool:
        res = pool.map(scan_illumina, [(f, per) for f in files])
    calls = [x for cl, _ in res for x in cl]
    n = max(1, len(calls))
    c = Counter(fam for fam, _, _ in calls)
    lens = sorted(ln for _, _, ln in calls)
    L = lens[-1] if lens else 0
    step = max(1, math.ceil(L / 100))
    curves = {}                          # family -> % of all reads with the adapter at or before each cycle bin
    for fam in ILLUMINA_ADAPTERS:
        at, run, cur = Counter(p // step for f, p, _ in calls if f == fam), 0, []
        for b in range(L // step + 1):
            run += at[b]
            cur.append(round(run / n, 5))
        curves[fam] = cur
    R = {"platform": "illumina", "mode": a.mode, "started": time.strftime("%Y-%m-%d %H:%M"), "n": len(calls),
         "files": len(files), "per_file": per, "found": ILLUMINA_FOUND, "readlen": lens[len(lens) // 2] if lens else 0, "maxlen": L, "step": step,
         "adapters": [[fam, lib, sq, c[fam] / n] for fam, (lib, sq, _) in ILLUMINA_ADAPTERS.items()] + [["none", "", "", c[None] / n]],
         "curves": curves}
    fam, m = max(((f, c[f]) for f in ILLUMINA_ADAPTERS), key=lambda kv: kv[1])
    if m >= ILLUMINA_FOUND * n:
        R["kit"] = fam
        R["also"] = [f for f in ILLUMINA_ADAPTERS if f != fam and c[f] >= ILLUMINA_FOUND * n]
        pos = sorted(p for f, p, _ in calls if f == fam)
        R["median_pos"] = pos[len(pos) // 2]
        R["verdict"] = f"the {fam} adapter is still on the reads with inserts shorter than the read (untrimmed)"
    else:
        R["verdict"] = "no known Illumina adapter in the reads: the inserts are longer than the reads, or the reads are trimmed"

    idx = sum((i for _, i in res), Counter())
    ni = max(1, sum(idx.values()))
    named = [(k, v) for k, v in idx.most_common() if k]
    if named:
        kind = lambda k: ("dual " if "+" in k else "single ") + "+".join(str(len(x)) for x in k.split("+"))
        kinds = Counter()
        for k, v in named:
            kinds[kind(k)] += v
        R["index"] = {"kinds": [[k, v / ni] for k, v in kinds.most_common() if v >= 0.01 * ni],
                      "top": [[k, v / ni] for k, v in named[:20]], "distinct": len(named), "none": idx[None] / ni}

    R["file_rows"] = []
    for f, (cl, ix) in zip(files, res):
        fc, fn = Counter(fam for fam, _, _ in cl), max(1, len(cl))
        row = {"file": label(f), "n": len(cl), "shares": {fam: fc[fam] / fn for fam in ILLUMINA_ADAPTERS}}
        top = [(k, v) for k, v in ix.most_common() if k]
        if top:
            row["idx"] = [top[0][0], top[0][1] / max(1, sum(ix.values()))]
        R["file_rows"].append(row)

    print(f"read context (illumina): {R['n']} reads from {len(files)} files; " +
          ", ".join(f"{k} {100 * f:.2f}%" for k, _, _, f in R["adapters"]) + f" -> {R.get('kit', 'none')}; index " +
          (", ".join(k for k, _ in R["index"]["kinds"]) if named else "none"), file=sys.stderr)
    write(a, R, ILLUMINA_BODY, ILLUMINA_JS)


HEAD = r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Read context</title>
<style>
:root{--bg:#fff;--panel:#fcfcfb;--panel2:#f1f0ec;--ink:#0b0b0b;--dim:#52514e;--muted:#898781;--grid:#e1e0d9;
  --ring:rgba(11,11,11,.10);--s1:#2a78d6;--good:#006300;--warn:#a86400;
  --q0:#f1f0ec;--q1:#cde2fb;--q2:#9ec5f4;--q3:#6da7ec;--q4:#3987e5;--q5:#256abf;
  --A:#1a9850;--C:#2a78d6;--G:#e6a100;--T:#d03b3b}
*{box-sizing:border-box}
html,body{margin:0;background:var(--bg);color:var(--ink);font:13px/1.45 system-ui,-apple-system,"Segoe UI",Roboto,sans-serif}
main{display:grid;gap:14px;padding:2px}
.panel{background:var(--panel);border:1px solid var(--ring);border-radius:10px;padding:12px 16px 14px;min-width:0}
h2{margin:-12px -16px 10px;padding:8px 16px;background:#e7eef8;border-bottom:1px solid var(--ring);border-radius:10px 10px 0 0;
  font-size:12px;letter-spacing:.06em;text-transform:uppercase;color:#1c3d66;font-weight:650}
h3{margin:14px 0 6px;font-size:12.5px;font-weight:600}
.verdict{font-size:14px;margin:0 0 6px}.verdict b{font-weight:650}.ok{color:var(--good)}.no{color:var(--warn)}
.two{display:grid;grid-template-columns:max-content minmax(0,1fr);gap:16px 40px;align-items:start}
@media (max-width:900px){.two{grid-template-columns:1fr}}
.scroll{overflow-x:auto}
table{border-collapse:collapse;font-variant-numeric:tabular-nums}
th,td{padding:4px 8px;border-bottom:1px solid var(--grid);text-align:right;white-space:nowrap;vertical-align:middle}
th{color:var(--dim);font-weight:600}th:first-child,td:first-child{text-align:left}td.l,th.l{text-align:left}
.seq{font-family:ui-monospace,Menlo,monospace;font-size:12px;letter-spacing:.02em}
.cell{display:inline-block;min-width:52px;border-radius:4px;padding:1px 6px;text-align:right}
tr.pick td{background:#e7eef8}
th.sort{cursor:pointer;user-select:none}th.sort:hover{color:var(--ink)}
.barsw{display:inline-block;margin-left:12px;font:12px system-ui;color:var(--dim);font-weight:normal}
.barsw button{font:12px system-ui;border:1px solid var(--ring);background:var(--panel);color:var(--dim);padding:1px 8px;cursor:pointer;margin-left:-1px}
.barsw button:first-child{margin-left:6px;border-radius:6px 0 0 6px}.barsw button:last-child{border-radius:0 6px 6px 0}
.barsw button.on{color:var(--ink);background:var(--panel2);font-weight:600}
.pbar{position:relative;height:14px;width:140px;background:var(--panel2);border:1px solid var(--grid)}
.pbar i{position:absolute;left:0;top:0;bottom:0;background:var(--s1)}.pbar.rel i{background:#d4831a}
.stack{display:flex;height:14px;width:220px;background:var(--panel2);border:1px solid var(--grid)}.stack i{display:block;height:100%}
.rv b{font-weight:650}.rv em{font-style:normal;color:var(--muted);margin-left:4px}
.note{color:var(--dim);font-size:12px;line-height:1.5;margin:6px 0}
.warn{color:var(--dim)}.warn::before{content:"! ";color:var(--warn);font-weight:700}
.legend{display:flex;flex-wrap:wrap;gap:4px 12px;color:var(--dim);font-size:11.5px;margin:4px 0}.legend span{white-space:nowrap}
.legend i{display:inline-block;width:10px;height:10px;border-radius:2px;margin-right:4px;vertical-align:-1px}
svg text{fill:var(--dim);font-size:10px}
</style></head><body>
"""

ONT_BODY = r"""<main>
  <section class="panel"><h2>Library kit at the read starts</h2><div id="verdict"></div>
    <div class="two long"><div id="kits"></div><div id="bcs"></div></div><div class="long" id="logo"></div></section>
  <section class="panel long"><h2>Per file</h2><div class="legend" id="fleg"></div><div class="scroll" id="files"></div></section>
</main>
"""

COMMON_JS = r"""const R = __DATA__;
const $ = id => document.getElementById(id);
const esc = s => String(s ?? "").replace(/[&<>"]/g, c => ({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;"}[c]));
const fmt = n => Number(n).toLocaleString("en-US");
const pf = f => (100 * f).toFixed(1) + "%";
const qv = ["--q0", "--q1", "--q2", "--q3", "--q4", "--q5"];
const heat = (f, lo = 0.0005) => { const i = f < lo ? 0 : Math.min(5, 1 + Math.floor(f * 5)); return `<span class="cell" style="background:var(${qv[i]});${i >= 4 ? "color:#fff" : ""}">${pf(f)}</span>`; };
const cseq = s => [...s].map(c => `<span style="color:var(--${"ACGT".includes(c) ? c : "muted"})">${esc(c)}</span>`).join("");
"""

ONT_JS = r"""const KC = {"SQK-NBD114": "#2a78d6", "SQK-LSK114": "#8fbcef", "SQK-MAB114": "#e08a00", "SQK-RBK114": "#8e57c4", "SQK-RAD114": "#cba6ea", "none": "#dcdbd5"};

// ---- verdict + kit table
$("verdict").innerHTML = (R.kit ? `<p class="verdict"><span class="ok">✓</span> read starts look like <b>${esc(R.kit)}</b>` +
                                  (R.mode === "short" ? "" : ` (${pf(R.kits.find(k => k[0] === R.kit)[2])} of reads)`) + `</p>`
                                : `<p class="verdict no">no known ONT kit at the read starts</p>`) +
  `<div class="note${R.kit ? "" : " warn"}">${esc(R.verdict)}</div>`;
if (R.mode === "short") document.querySelectorAll(".long").forEach(e => e.style.display = "none");
$("kits").innerHTML = `<h3>kit adapters in the first ${R.window} nt</h3><table><tr><th class="l">kit</th><th class="l">read start</th><th>reads</th></tr>` +
  R.kits.map(([k, lay, f]) => `<tr class="${k === R.kit ? "pick" : ""}"><td class="l"><b>${esc(k)}</b></td><td class="l">${esc(lay || "no known adapter")}</td><td>${heat(f)}</td></tr>`).join("") +
  `</table><div class="note">${fmt(R.n)} reads from ${fmt(R.files)} files (up to ${fmt(R.per_file)} per file); known adapters within ~15% edits</div>`;

// ---- base frequencies of the kit reads aligned at the insert, kit layout as bands behind the letters
function logo(e, lay, bcs) {
  const U = e.logo.length, cw = 12, lw = 44, H = 48, Wd = lw + U * cw + 46, top = 24, base = top + H;
  const xi = p => lw + (U - p) * cw;      // left edge of position -p (p = 1 is the base right before the insert)
  let h = `<svg viewBox="0 0 ${Wd} ${base + 58}" width="${Wd}" style="max-width:100%;height:auto" role="img" aria-label="base frequencies before the insert">`;
  const bc = bcs.length && bcs[0][1] >= 0.8 ? bcs[0] : null;           // one barcode: show its sequence
  const bcLab = bc ? ` (${bc[0]})` : bcs.length > 1 ? ` (mixed: ${bcs.slice(0, 3).map(b => b[0]).join(", ")}${bcs.length > 3 ? ", …" : ""})` : "";
  let p = 0, alt = 0, bcSeq = "";
  [...lay].reverse().forEach(([lab, sq, n]) => {
    const p0 = p + 1, p1 = p + n; p = p1;
    bcSeq = (sq || (bc ? bc[2] : "·".repeat(n))) + bcSeq;
    if (p0 > U) return;
    const x0 = xi(Math.min(p1, U)), x1 = xi(p0) + cw, isBc = !sq;
    const fill = isBc ? "#fdf0d5" : alt++ % 2 ? "#e2e9f2" : "#eef2f7";
    h += `<rect x="${x0 + 1}" y="${top - 16}" width="${x1 - x0 - 2}" height="${H + 16}" fill="${fill}" rx="3"/>`;
    h += `<text x="${(x0 + x1) / 2}" y="${top - 5}" text-anchor="middle" style="fill:${isBc ? "#8a5a00" : "var(--dim)"};font-weight:${isBc ? 650 : 500}">${esc(lab)}${isBc ? esc(bcLab) : ""}${p1 > U ? " …" : ""}</text>`;
  });
  for (const f of [0, .5, 1]) h += `<line x1="${lw}" x2="${lw + U * cw}" y1="${base - f * H}" y2="${base - f * H}" stroke="var(--grid)"/><text x="${lw - 4}" y="${base - f * H + 3}" text-anchor="end">${f * 100}%</text>`;
  e.logo.forEach((fr, i) => { let y = base;
    ["A", "C", "G", "T"].map((b, k) => [b, fr[k]]).sort((a, b) => a[1] - b[1]).forEach(([b, f]) => { const hh = f * H; if (hh < 1) return;
      h += `<text transform="translate(${lw + i * cw + 0.5} ${y.toFixed(2)}) scale(${(cw - 2) / 7.2} ${(hh / 9.4).toFixed(3)})" style="font:400 13px ui-monospace,Menlo,monospace;fill:var(--${b})">${b}<title>-${U - i}: ${b} ${pf(f)}</title></text>`; y -= hh; }); });
  h += `<text x="${lw + U * cw + 4}" y="${base - 2}" style="fill:var(--ink);font-weight:600">insert</text>`;
  for (let i = 0; i < U; i += 10) h += `<text x="${lw + i * cw + cw / 2}" y="${base + 12}" text-anchor="middle">-${U - i}</text>`;
  const row = (y, lab, s) => { let t = `<text x="${lw - 4}" y="${y}" text-anchor="end">${lab}</text>`;
    [...s].forEach((c, i) => { const x = lw + (U - s.length + i) * cw + cw / 2; if (x < lw) return;
      t += `<text x="${x}" y="${y}" text-anchor="middle" style="font:600 11px ui-monospace,Menlo,monospace;fill:var(--${"ACGT".includes(c) ? c : "muted"})">${esc(c)}</text>`; });
    return t; };
  h += row(base + 30, "cons.", e.consensus);
  h += row(base + 48, "kit", bcSeq);
  return h + "</svg>";
}
if (R.align) $("logo").innerHTML = `<h3>${esc(R.kit)} reads aligned at the insert</h3><div class="legend">${["A", "C", "G", "T"].map(b => `<span><i style="background:var(--${b})"></i>${b}</span>`).join("")}` +
  `<span style="white-space:normal">letter height = fraction of the ${fmt(R.align.reads)} ${esc(R.kit)} reads with that base</span></div>` + logo(R.align, R.align.layout, R.barcodes || []);

// ---- ONT barcodes, sortable by barcode or count; top BC_TOP rows, all found, or the kit's whole set
const BT = {}, BC_TOP = 6;
function bcTable(id, rows, checked, note, set) {
  const got = new Set(rows.map(r => r[0]));
  BT[id] = {rows, checked, note, key: "n", dir: -1, view: "top", rel: false,
            set: rows.concat((set || []).filter(([b]) => !got.has(b)).map(([b, sq]) => [b, 0, sq]))};
  drawBc(id);
}
function viewBc(id, view) { BT[id].view = view; drawBc(id); }
function barBc(id, rel) { BT[id].rel = rel; drawBc(id); }
function sortBc(id, key) { const t = BT[id]; t.dir = t.key === key ? -t.dir : key === "b" ? 1 : -1; t.key = key; drawBc(id); }
function drawBc(id) {
  const t = BT[id], top = Math.max(1e-9, ...t.rows.map(r => r[1]));
  const rows = [...(t.view === "set" ? t.set : t.rows)].sort((a, b) => t.dir * (t.key === "b" ? a[0].localeCompare(b[0], "en", {numeric: true}) : a[1] - b[1]));
  const th = (key, label, attr) => `<th class="sort${key === "n" ? " l" : ""}"${attr} onclick="sortBc('${id}','${key}')">${label}${t.key === key ? (t.dir > 0 ? " ▲" : " ▼") : ""}</th>`;
  const vw = (v, label) => `<button class="${t.view === v ? "on" : ""}" onclick="viewBc('${id}','${v}')">${label}</button>`;
  const sw = (rel, label) => `<button class="${t.rel === rel ? "on" : ""}" onclick="barBc('${id}',${rel})">${label}</button>`;
  $(id).innerHTML = `<h3>ONT barcodes · ${rows.length}` +
    `<span class="barsw">show ${vw("top", "top " + BC_TOP)}${vw("found", "found (" + t.rows.length + ")")}${vw("set", "whole set (" + t.set.length + ")")}</span>` +
    `<span class="barsw">bar ${sw(false, "share of all")}${sw(true, "share of top")}</span></h3><table><tr>${th("b", "barcode", "")}<th class="l">sequence</th>${th("n", "reads", ' colspan="2"')}</tr>` +
    rows.slice(0, t.view === "top" ? BC_TOP : rows.length).map(([b, f, sq]) => `<tr><td>${esc(b)}</td><td class="l seq">${cseq(sq || "")}</td>` +
      `<td class="l"><div class="pbar${t.rel ? " rel" : ""}"><i style="width:${(100 * (t.rel ? f / top : f)).toFixed(2)}%"></i></div></td>` +
      `<td class="l rv"><b>${pf(f)}</b> <em>${fmt(Math.round(f * t.checked))}</em></td></tr>`).join("") +
    `</table><div class="note">${t.note}</div>`;
}
if (R.barcodes)
  bcTable("bcs", R.barcodes, R.checked, `of ${fmt(R.checked)} ${esc(R.kit)} reads checked (a few per file); ${pf(R.unclear)} unclear`, R.set);

// ---- per file: kit shares as a stacked bar, main kit, barcodes
const order = R.kits.map(k => k[0]);
$("fleg").innerHTML = order.map(k => `<span><i style="background:${KC[k]}"></i>${esc(k)}</span>`).join("");
$("files").innerHTML = `<table><tr><th class="l">file</th><th>reads</th><th class="l">kit at the read start</th><th class="l">main</th>${R.barcodes ? `<th class="l">ONT barcode</th>` : ""}</tr>` +
  R.file_rows.map(r => { const main = Object.entries(r.shares).sort((a, b) => b[1] - a[1])[0] || ["none", 0];
    return `<tr><td class="l">${esc(r.file)}</td><td>${fmt(r.n)}</td>` +
      `<td class="l"><div class="stack" title="${order.filter(k => r.shares[k]).map(k => k + " " + pf(r.shares[k])).join(", ")}">${order.map(k => r.shares[k] ? `<i style="width:${(100 * r.shares[k]).toFixed(2)}%;background:${KC[k]}"></i>` : "").join("")}</div></td>` +
      `<td class="l rv"><b>${esc(main[0])}</b> <em>${pf(main[1])}</em></td>` +
      (R.barcodes ? `<td class="l rv" title="${fmt(r.bc_n || 0)} reads checked">${r.bc && r.bc.length ? r.bc.map(([b, f]) => `<b>${esc(b)}</b> <em>${pf(f)}</em>`).join(", ") : `<em>–</em>`}</td>` : "") + `</tr>`; }).join("") +
  `</table>` + (R.barcodes ? `<div class="note">ONT barcode: up to ${fmt(Math.max(...R.file_rows.map(r => r.bc_n || 0)))} ${esc(R.kit)} reads checked per file (hover for the count)</div>` : "");

"""

ILLUMINA_BODY = r"""<main>
  <section class="panel"><h2>Adapter read-through</h2><div id="verdict"></div>
    <div class="long" id="ads"></div><div class="long" id="curve"></div></section>
  <section class="panel long"><h2>Indexes in the read headers</h2><div id="idx"></div></section>
  <section class="panel long"><h2>Per file</h2><div class="scroll" id="files"></div></section>
</main>
"""

ILLUMINA_JS = r"""const AC = {"TruSeq": "#2a78d6", "Nextera": "#e08a00", "small RNA": "#8e57c4", "none": "#dcdbd5"};
const fams = R.adapters.filter(r => r[0] !== "none").map(r => r[0]);

// ---- verdict + adapter table
$("verdict").innerHTML = (R.kit ? `<p class="verdict"><span class="ok">✓</span> adapter read-through: <b>${esc(R.kit)}</b>${R.also.length ? ` (also ${R.also.map(esc).join(", ")})` : ""}` +
                                  (R.mode === "short" ? "" : ` (${pf(R.adapters.find(r => r[0] === R.kit)[3])} of reads, median start at cycle ${fmt(R.median_pos + 1)})`) + `</p>`
                                : `<p class="verdict">no adapter read-through</p>`) +
  `<div class="note">${esc(R.verdict)}</div>`;
if (R.mode === "short") document.querySelectorAll(".long").forEach(e => e.style.display = "none");
$("ads").innerHTML = `<h3>adapters anywhere in the read</h3><table><tr><th class="l">adapter</th><th class="l">library</th><th class="l">sequence</th><th>reads</th></tr>` +
  R.adapters.map(([k, lib, sq, f]) => `<tr class="${k === R.kit ? "pick" : ""}"><td class="l"><b>${esc(k)}</b></td><td class="l">${esc(lib || "no known adapter")}</td>` +
    `<td class="l seq">${cseq(sq)}</td><td>${heat(f, k === "none" ? 0.0005 : R.found)}</td></tr>`).join("") +
  `</table><div class="note">${fmt(R.n)} reads from ${fmt(R.files)} files (up to ${fmt(R.per_file)} per file), read length ${fmt(R.readlen)} (median); adapter start within 1-2 edits</div>`;

// ---- cumulative adapter content by cycle (as in FastQC): % of reads with the adapter starting at or before a cycle
{ const W = 560, H = 160, l = 44, b = 22, t = 8, r = 120, top = Math.max(0.01, ...fams.map(f => R.curves[f].at(-1) || 0));
  const ymax = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1].find(v => v >= top) || 1, nb = R.curves[fams[0]].length;
  const x = i => l + (nb > 1 ? i / (nb - 1) : 0) * (W - l - r), y = v => t + (1 - v / ymax) * (H - t - b);
  let h = `<svg viewBox="0 0 ${W} ${H}" width="${W}" style="max-width:100%;height:auto" role="img" aria-label="adapter content by cycle">`;
  for (const f of [0, .5, 1]) h += `<line x1="${l}" x2="${W - r}" y1="${y(f * ymax)}" y2="${y(f * ymax)}" stroke="var(--grid)"/><text x="${l - 4}" y="${y(f * ymax) + 3}" text-anchor="end">${pf(f * ymax)}</text>`;
  for (const f of [0, .25, .5, .75, 1]) { const i = Math.round(f * (nb - 1)); h += `<text x="${x(i)}" y="${H - 6}" text-anchor="middle">${fmt(i * R.step + 1)}</text>`; }
  fams.forEach(f => { const c = R.curves[f], last = c.at(-1) || 0;
    h += `<polyline fill="none" stroke="${AC[f]}" stroke-width="2" points="${c.map((v, i) => `${x(i).toFixed(1)},${y(v).toFixed(1)}`).join(" ")}"><title>${esc(f)}</title></polyline>`;
    h += `<text x="${W - r + 4}" y="${y(last) + 3}" style="fill:${AC[f]};font-weight:600">${esc(f)} ${pf(last)}</text>`; });
  $("curve").innerHTML = `<h3>adapter content by cycle</h3>${h}</svg><div class="note">% of reads with the adapter starting at or before the cycle; ` +
    `an early rise means many inserts shorter than the read</div>`; }

// ---- indexes
const I = R.index, iseq = s => s.split("+").map(cseq).join(`<span style="color:var(--muted)"> + </span>`);
$("idx").innerHTML = !I ? `<div class="note">no index in the read headers (CASAVA 1.8+ format, e.g. <span class="seq">1:N:0:ATCACG+TTAGGC</span>)</div>` :
  `<p class="verdict">${I.kinds.map(([k, f]) => `<b>${esc(k)}</b> nt` + (I.kinds.length > 1 ? ` <span style="color:var(--muted)">${pf(f)}</span>` : "")).join(", ")} · ${fmt(I.distinct)} distinct indexes</p>` +
  `<table><tr><th class="l">index</th><th class="l" colspan="2">reads</th></tr>` +
  I.top.filter(([, f]) => f >= 0.005).slice(0, 10).map(([k, f]) => `<tr><td class="l seq">${iseq(k)}</td><td class="l"><div class="pbar"><i style="width:${(100 * f / I.top[0][1]).toFixed(2)}%"></i></div></td><td class="l rv"><b>${pf(f)}</b></td></tr>`).join("") +
  `</table><div class="note">indexes in ≥ 0.5% of reads (up to 10), as read by the sequencer; bar = share of the top index` +
  `${I.none > 0 ? `; ${pf(I.none)} of reads without an index in the header` : ""}</div>`;

// ---- per file
$("files").innerHTML = `<table><tr><th class="l">file</th><th>reads</th>${fams.map(f => `<th>${esc(f)}</th>`).join("")}${I ? `<th class="l">index</th>` : ""}</tr>` +
  R.file_rows.map(r => `<tr><td class="l">${esc(r.file)}</td><td>${fmt(r.n)}</td>${fams.map(f => `<td>${heat(r.shares[f] || 0, R.found)}</td>`).join("")}` +
    (I ? `<td class="l rv">${r.idx ? `<span class="seq">${iseq(r.idx[0])}</span> <em>${pf(r.idx[1])}</em>` : `<em>–</em>`}</td>` : "") + `</tr>`).join("") +
  `</table><div class="note">% of reads with the adapter (shaded from ${pf(R.found)}, below that it is noise); index = most common index in the file's read headers</div>`;
"""

FIT_JS = r"""// grow the iframe (faster-report shows this page in one) to the content
function fit() { try { const f = window.frameElement; if (f) f.style.height = Math.ceil(document.querySelector("main").getBoundingClientRect().bottom + 4) + "px"; } catch (e) {} }
addEventListener("load", fit); addEventListener("resize", fit); fit();
"""


def page(body, js):
    return HEAD + body + "<script>\n" + COMMON_JS + js + FIT_JS + "</script></body></html>\n"

if __name__ == "__main__":
    main()
