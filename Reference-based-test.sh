#!/usr/bin/env bash
# CORONAS: GDC2 reference-based compression vs published reference-free OptimJV3.
# Run on Ubuntu: bash run_CORONAS_GDC2_all.sh
# The GDC2 archive contains several files; ALL components are accounted for.
set -Eeuo pipefail
export LC_ALL=C

BASE="${WORKDIR:-$(pwd)/coronas_GDC2_benchmark}"
mkdir -p "$BASE"
BASE="$(cd "$BASE" && pwd)"
DATA="$BASE/data"; REPO="$BASE/repos/GDC2"; RESULTS="$BASE/results"; LOGS="$BASE/logs"
mkdir -p "$DATA" "$BASE/repos" "$RESULTS" "$LOGS"
CORONAS_URL="https://raw.githubusercontent.com/cobilab/OptJV3/main/coronas.fasta.zst"
REF_ACCESSION=NC_045512.2
REF_URL="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${REF_ACCESSION}&rettype=fasta&retmode=text"
EXPECTED_SYMBOLS="${EXPECTED_SYMBOLS:-58121050}"
OPT_BYTES=85050; JARVIS_BYTES=112387
THREADS="${THREADS:-6}"; MEMORY_GB="${MEMORY_GB:-4}"; TIMEOUT_S="${TIMEOUT_S:-7200}"
# Internal GDC2 allocation limit (MiB); must be lower than the process virtual-memory cap.
GDC2_MM="${GDC2_MM:-3072}"

info(){ printf '\n[%s] %s\n' "$(date +%H:%M:%S)" "$*"; }
die(){ printf 'ERROR: %s\n' "$*" >&2; exit 1; }
trap 'printf "ERROR at line %s (command exit %s). Check %s\n" "$LINENO" "$?" "$LOGS" >&2' ERR

need_packages(){
  # Explicitly check zlib1g-dev, even when git/make/g++ are already installed.
  # The GDC2 linker uses zlib's gzopen/gzwrite/gzread symbols.
  local needed=(git build-essential make curl zstd python3 time coreutils ca-certificates zlib1g-dev)
  local missing=() pkg
  for pkg in "${needed[@]}"; do
    if ! dpkg-query -W -f='${Status}' "$pkg" 2>/dev/null | grep -q 'install ok installed'; then
      missing+=("$pkg")
    fi
  done
  if ((${#missing[@]})); then
    command -v apt-get >/dev/null || die "Missing packages: ${missing[*]} (requires Ubuntu/Debian)"
    local elevate=()
    if (( EUID != 0 )); then
      command -v sudo >/dev/null || die "sudo needed to install packages: ${missing[*]}"
      elevate=(sudo)
    fi
    info "Installing dependencies (including zlib1g-dev): ${missing[*]}"
    "${elevate[@]}" apt-get update
    "${elevate[@]}" apt-get install -y "${missing[@]}"
  else
    info 'Dependencies already installed, including zlib1g-dev'
  fi
  [[ -e /usr/include/zlib.h ]] || die 'zlib.h not found despite package installation'
}

get_file(){
  local url="$1" dest="$2" tmp="${2}.part"
  if [[ -s "$dest" ]]; then info "Cached: $dest"; return; fi
  rm -f "$tmp"
  curl --fail --location --silent --show-error --retry 6 --retry-all-errors --retry-delay 3 \
       --connect-timeout 30 --max-time 900 -o "$tmp" "$url"
  [[ -s "$tmp" ]] || die "Empty download: $url"
  mv "$tmp" "$dest"
}

need_packages
info '1/5 Downloading the exact CORONAS collection used in the manuscript'
if [[ -n "${CORONAS_FASTA:-}" ]]; then
  [[ -s "$CORONAS_FASTA" ]] || die "CORONAS_FASTA not found: $CORONAS_FASTA"
  # Source path can be identical to destination if rerunning.
  [[ "$(realpath "$CORONAS_FASTA")" == "$(realpath -m "$DATA/coronas.fasta")" ]] || cp "$CORONAS_FASTA" "$DATA/coronas.fasta"
else
  get_file "$CORONAS_URL" "$DATA/coronas.fasta.zst"
  zstd --quiet --test "$DATA/coronas.fasta.zst" || die 'CORONAS .zst archive failed verification'
  [[ -s "$DATA/coronas.fasta" ]] || zstd -dc "$DATA/coronas.fasta.zst" > "$DATA/coronas.fasta"
fi
info "2/5 Downloading the SARS-CoV-2 reference $REF_ACCESSION"
if [[ -n "${REFERENCE_FASTA:-}" ]]; then
  [[ -s "$REFERENCE_FASTA" ]] || die "REFERENCE_FASTA not found: $REFERENCE_FASTA"
  [[ "$(realpath "$REFERENCE_FASTA")" == "$(realpath -m "$DATA/${REF_ACCESSION}.fasta")" ]] || cp "$REFERENCE_FASTA" "$DATA/${REF_ACCESSION}.fasta"
else
  get_file "$REF_URL" "$DATA/${REF_ACCESSION}.fasta"
fi

info '3/5 Preparing one normalized FASTA file per CORONAS genome (GDC2 requirement)'
# Importantly, GDC2 CANNOT use the single concatenated .seq input: it needs the
# genome boundaries. The order and ACGT filtering match OptJV3 preprocessing.
python3 - "$DATA/coronas.fasta" "$DATA/${REF_ACCESSION}.fasta" "$DATA/processed" "$EXPECTED_SYMBOLS" "$RESULTS/dataset.json" <<'PY'
import hashlib, json, re, sys
from pathlib import Path
original,reference,output,expected,manifest=(Path(sys.argv[1]),Path(sys.argv[2]),Path(sys.argv[3]),int(sys.argv[4]),Path(sys.argv[5]))
output.mkdir(parents=True,exist_ok=True)

def records(file):
    header=None; parts=[]
    with file.open('rt',encoding='ascii',errors='replace') as f:
        for line in f:
            if line.startswith('>'):
                if header is not None: yield header,''.join(parts)
                header=line[1:].strip(); parts=[]
            elif header is not None:
                parts.append(line.strip().upper())
    if header is not None: yield header,''.join(parts)

def normalize(seq): return re.sub('[^ACGT]','',seq.upper())
ref=list(records(reference))
if len(ref)!=1: raise SystemExit('Expected exactly one reference FASTA record')
ref_seq=normalize(ref[0][1]); nref=len(ref_seq)
if not 28000<nref<31000: raise SystemExit(f'Unexpected reference size: {nref}')
(output/'reference.fa').write_text('>NC_045512.2\n'+ref_seq+'\n',encoding='ascii')
files=[]; n=0; original_count=0; symbols=0; lost=0; sha=hashlib.sha256(); sizes=[]
for header,seq in records(original):
    original_count+=1
    cleaned=normalize(seq)
    lost+=len(seq)-len(cleaned)
    if not cleaned: raise SystemExit(f'Empty normalized genome: record {original_count}')
    n+=1; symbols+=len(cleaned); sha.update(cleaned.encode('ascii'))
    name=f'genome_{n:05d}.fa'; files.append(name); sizes.append(len(cleaned))
    (output/name).write_text(f'>genome_{n:05d}\n{cleaned}\n',encoding='ascii')
if n<2: raise SystemExit(f'GDC2 needs individual genomes; found only {n} FASTA records. Cannot infer boundaries from concatenated data.')
if symbols!=expected: raise SystemExit(f'CORONAS ACGT length {symbols} differs from manuscript value {expected}; refusing invalid comparison')
(output/'files.txt').write_text('reference.fa\n'+'\n'.join(files)+'\n',encoding='ascii')
meta={'dataset':'CORONAS','genomes':n,'source_records':original_count,'acgt_symbols':symbols,'removed_non_acgt':lost,'reference_accession':'NC_045512.2','reference_acgt_symbols':nref,'sha256_concatenated_acgt':sha.hexdigest(),'reference_sha256':hashlib.sha256(ref_seq.encode()).hexdigest(),'min_genome_length':min(sizes),'max_genome_length':max(sizes)}
manifest.write_text(json.dumps(meta,indent=2)+'\n')
print(f'CORONAS: {n} genomes, {symbols:,} ACGT symbols; reference: {nref:,} ACGT; removed ambiguous symbols: {lost:,}')
print('Concatenated ACGT SHA256:',sha.hexdigest())
PY

info '4/5 Installing/building GDC2 (refresh-bio/GDC2)'
if [[ -n "${GDC2_BIN:-}" ]]; then
  GDC2="$(realpath "$GDC2_BIN")"
  [[ -x "$GDC2" ]] || die "GDC2_BIN is not executable: $GDC2"
  GDC_COMMIT=provided_binary
else
  if [[ ! -d "$REPO/.git" ]]; then
    git clone --depth 1 https://github.com/refresh-bio/GDC2.git "$REPO"
  fi
  GDC_COMMIT="$(git -C "$REPO" rev-parse HEAD)"
  MAKEFILE="$REPO/gdc_2/Gdc2/makefile.linux"
  [[ -f "$MAKEFILE" ]] || die "Missing GDC2 makefile: $MAKEFILE"
  # Fix GDC2 upstream's zlib link order. GNU ld processes -lz when seen;
  # if it appears *before* the .o files, gzopen/gzread/etc. remain undefined.
  # This patch is idempotent: safe to rerun after an interrupted benchmark.
  python3 - "$MAKEFILE" <<'PATCH_MAKEFILE'
import re, sys
from pathlib import Path
p=Path(sys.argv[1]); original=p.read_text()
lines=original.splitlines(keepends=True)
found_flags=False; found_link=False
for i,line in enumerate(lines):
    if re.match(r'^CLINK\s*=',line):
        lines[i]=re.sub(r'(?<!\S)-lz(?=\s|$)','',line).rstrip()+'\n'
        found_flags=True
    elif '$(CC) $(CLINK) -o gdc2' in line:
        if '../libs/libaelf64.a' not in line:
            raise SystemExit('Unexpected GDC2 link command; cannot patch safely')
        line=re.sub(r'(?:\s+-(?:lpthread|lm|lz))+\s*$', '', line).rstrip()
        # Keep pthread, math and zlib AFTER all object files and static libs.
        lines[i]=line+' -lpthread -lm -lz\n'
        found_link=True
if not (found_flags and found_link):
    raise SystemExit('Unrecognized GDC2 makefile layout; refusing unsafe build')
new=''.join(lines)
if new!=original:
    p.write_text(new)
    print('Patched makefile.linux: linked -lpthread -lm -lz AFTER object files')
else:
    print('GDC2 makefile.linux already patched')
PATCH_MAKEFILE
  if [[ ! -x "$REPO/gdc_2/Gdc2/gdc2" ]]; then
    info "Building GDC2; build log: $LOGS/build_GDC2.log"
    if ! make -C "$REPO/gdc_2/Gdc2" -f makefile.linux -j4 > "$LOGS/build_GDC2.log" 2>&1; then
      tail -n 55 "$LOGS/build_GDC2.log" >&2
      die "GDC2 build failed; consult $LOGS/build_GDC2.log"
    fi
  else
    info 'Using previously built GDC2 binary'
  fi
  GDC2="$REPO/gdc_2/Gdc2/gdc2"
fi
[[ -x "$GDC2" ]] || die "Cannot find GDC2 executable: $GDC2"
[[ "$THREADS" =~ ^[0-9]+$ ]] && ((THREADS>=2)) || die 'THREADS must be >= 2'
MEM_KB="$(python3 - "$MEMORY_GB" <<'PY'
import sys
print(int(float(sys.argv[1])*1024*1024))
PY
)"
[[ "$GDC2_MM" =~ ^[0-9]+$ ]] || die 'GDC2_MM must be positive integer'

# The first item is the reference and is encoded in *.gdc2_ref; later items
# are the CORONAS genomes. The archive has MULTIPLE component files.
ARCHIVE="$RESULTS/coronas_gdc2"
# By default, reuse an existing complete archive. This allows recovery from a
# failed verification without wasting time recompressing all 2,000 genomes.
if [[ "${FORCE_RECOMPRESS:-0}" != "1" && -s "${ARCHIVE}.gdc2_ref" && -s "${ARCHIVE}.gdc2_desc" && -s "${ARCHIVE}.gdc2_rc" && -s "$RESULTS/gdc2_metrics.txt" ]]; then
  info "5/5 Reusing existing GDC2 archive (set FORCE_RECOMPRESS=1 to recompress)"
else
  rm -f "${ARCHIVE}".gdc2_*
  info "5/5 Compressing ${EXPECTED_SYMBOLS} symbols with GDC2 (-t${THREADS}, -mm${GDC2_MM})"
  info "Limit: ${MEMORY_GB} GiB virtual memory per GDC2 process; timeout: ${TIMEOUT_S}s"
  set +e
  (cd "$DATA/processed" && /usr/bin/time -f 'WALL_S=%e\nMAXRSS_KB=%M' -o "$RESULTS/gdc2_metrics.txt" \
    timeout --signal=TERM --kill-after=10s "$TIMEOUT_S" bash -c \
    'ulimit -v "$1"; shift; exec "$@"' bash "$MEM_KB" \
    "$GDC2" c "-t${THREADS}" "-mm${GDC2_MM}" "$ARCHIVE" '@files.txt' \
    > "$LOGS/GDC2_compression.log" 2>&1)
  status=$?
  set -e
  if ((status!=0)); then
    tail -n 45 "$LOGS/GDC2_compression.log" >&2
    die "GDC2 exited ${status}; see $LOGS/GDC2_compression.log"
  fi
fi
for suffix in ref desc rc; do
  [[ -s "${ARCHIVE}.gdc2_${suffix}" ]] || die "Missing archive component ${ARCHIVE}.gdc2_${suffix}"
done

# GDC2 writes decoded genomes as <original filename>.ori (not <original filename>).
# Empty the dedicated verify directory before each verification to rule out stale outputs.
VERIFY="$RESULTS/verify"
mkdir -p "$VERIFY"
find "$VERIFY" -maxdepth 1 -type f -name '*.ori' -delete
info "Checking GDC2 round-trip for all genomes (verify log: $LOGS/GDC2_verify.log)"
set +e
(cd "$VERIFY" && timeout --signal=TERM --kill-after=10s "$TIMEOUT_S" "$GDC2" d "$ARCHIVE" \
  > "$LOGS/GDC2_verify.log" 2>&1)
verify_status=$?
set -e
if ((verify_status!=0)); then
  tail -n 45 "$LOGS/GDC2_verify.log" >&2
  die "GDC2 decompression failed (exit ${verify_status}); see $LOGS/GDC2_verify.log"
fi
python3 - "$DATA/processed" "$VERIFY" <<'PY_VERIFY'
from pathlib import Path
import hashlib, sys
src, decoded = map(Path, sys.argv[1:])
expected = [line.strip() for line in (src / 'files.txt').read_text().splitlines() if line.strip()]

def normalized_sequence(path):
    if not path.is_file():
        raise SystemExit(f'Missing decompressed file: {path} (GDC2 should append .ori)')
    h = hashlib.sha256()
    count = 0
    with path.open('rb') as stream:
        for line in stream:
            if line.startswith(b'>'):
                continue
            seq = line.strip().upper()
            if any(b not in b'ACGT' for b in seq):
                raise SystemExit(f'Unexpected non-ACGT symbols in {path}')
            h.update(seq)
            count += len(seq)
    return count, h.hexdigest()

for i, name in enumerate(expected):
    original = src / name
    recovered = decoded / (name + '.ori')
    if normalized_sequence(original) != normalized_sequence(recovered):
        raise SystemExit(f'Decompression mismatch for {name}')
    if i and i % 250 == 0:
        print(f'[verification] {i}/{len(expected)-1} genomes checked', flush=True)
print(f'PASS: all {len(expected)-1} CORONAS genomes and the reference recovered losslessly')
PY_VERIFY

# Size includes the compressed reference, descriptors, and all GDC2 components.
python3 - "$ARCHIVE" "$RESULTS/gdc2_metrics.txt" "$RESULTS/dataset.json" "$RESULTS/summary.tsv" "$RESULTS/summary.md" "$OPT_BYTES" "$JARVIS_BYTES" "$GDC_COMMIT" "$THREADS" "$MEMORY_GB" "$GDC2_MM" <<'PY'
import json,sys
from pathlib import Path
archive,metrics,manifest,tsv,md=map(Path,sys.argv[1:6]);opt,jarvis=map(int,sys.argv[6:8]);commit=sys.argv[8];threads=sys.argv[9];mem=sys.argv[10];mm=sys.argv[11]
data=json.loads(manifest.read_text()); n=data['acgt_symbols']
components=sorted(archive.parent.glob(archive.name+'.gdc2_*'))
if not components: raise SystemExit('No GDC2 archive files found')
full=sum(p.stat().st_size for p in components)
ref=archive.parent/(archive.name+'.gdc2_ref')
if not ref.exists(): raise SystemExit('No GDC2 reference archive found')
ref_bytes=ref.stat().st_size
conditional=full-ref_bytes
mt={}
for line in metrics.read_text().splitlines():
    if '=' in line:
        k,v=line.strip().split('=',1);mt[k]=v
wall=mt.get('WALL_S','');rss=round(int(mt.get('MAXRSS_KB','0'))/1048576,5)
rows=[
 ['GDC2_full_archive','measured_verified',full,f'{8*full/n:.8f}',wall,rss,'yes', f'{100*(full-opt)/opt:+.3f}'],
 ['GDC2_excluding_reference_component','conditional_not_standalone',conditional,f'{8*conditional/n:.8f}',wall,rss,'no', f'{100*(conditional-opt)/opt:+.3f}'],
 ['OptimJV3_manuscript','published_not_rerun',opt,f'{8*opt/n:.8f}','','','no','0.000'],
 ['JARVIS3_level33_manuscript','published_not_rerun',jarvis,f'{8*jarvis/n:.8f}','','','no',f'{100*(jarvis-opt)/opt:+.3f}']]
heads=['method','status','compressed_bytes','bps','compression_time_s','peak_rss_gb','reference_included','size_gap_vs_OptimJV3_pct']
tsv.write_text('\t'.join(heads)+'\n'+'\n'.join('\t'.join(map(str,r)) for r in rows)+'\n')
with md.open('w') as f:
    f.write('| '+' | '.join(heads)+' |\n| '+' | '.join(['---']*len(heads))+' |\n')
    for r in rows: f.write('| '+' | '.join(map(str,r))+' |\n')
    f.write('\n*The GDC2 conditional figure excludes the `.gdc2_ref` archive component and is not independently decodable; use the full-archive number for self-contained storage. Published JARVIS3/OptimJV3 values were not rerun.*\n')
json_out=tsv.with_name('environment.json')
json_out.write_text(json.dumps({'gdc2_git_commit':commit,'threads':threads,'process_virtual_memory_gib':mem,'gdc2_internal_memory_mib':mm,'source_dataset':data,'archive_components':{p.name:p.stat().st_size for p in components},'total_archive_bytes':full,'reference_component_bytes':ref_bytes,'conditional_target_bytes':conditional,'comparison_reference_free_bytes':opt,'compression_wall_s':wall,'compression_peak_rss_gib':rss},indent=2)+'\n')
print(md.read_text())
print('Archive components:',[(p.name,p.stat().st_size) for p in components])
print('Results:',tsv,md)
PY
info "Benchmark complete: $RESULTS/summary.tsv and $RESULTS/summary.md"
info 'IMPORTANT: GDC2 full-archive includes its reference; conditional figure excludes the compressed reference component.'
info 'IMPORTANT: 85,050 B (OptJV3) and 112,387 B (JARVIS3) are previously published, NOT re-measured.'
