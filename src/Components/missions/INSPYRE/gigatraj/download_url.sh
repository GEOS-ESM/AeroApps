#!/usr/bin/env bash

set -Eeuo pipefail

usage() {
  cat <<'EOF'
Usage:
  download_url_glob.sh [options] URL GLOB

Download files linked by an HTTP(S) directory listing whose decoded filenames
match a shell glob.

Options:
  -o, --output DIR       Destination directory (default: .)
  -j, --jobs N           Concurrent downloads (default: 4)
      --skip N           Skip the first N files in the ordered match list
      --numFiles N       Download at most N files after skipping
      --tool TOOL        auto, curl, or wget (default: auto)
  -n, --dry-run          Print matches without downloading
  -h, --help             Show this help

Example:
  ./download_url_glob.sh -j 6 --skip 4 --numFiles 8 -o ./asm \
    'https://portal.nccs.nasa.gov/datashare/gmao/geos-fp/forecast/Y2026/M08/D28/H00/' \
    '*_asm_Nv.*.nc4'
EOF
}

output_dir=.
jobs=4
skip=0
num_files=-1
tool=auto
dry_run=0

while (($#)); do
  case "$1" in
    -o|--output) output_dir=${2:?"missing value for $1"}; shift 2 ;;
    -j|--jobs) jobs=${2:?"missing value for $1"}; shift 2 ;;
    --skip) skip=${2:?"missing value for $1"}; shift 2 ;;
    --numFiles|--num-files) num_files=${2:?"missing value for $1"}; shift 2 ;;
    --tool) tool=${2:?"missing value for $1"}; shift 2 ;;
    -n|--dry-run) dry_run=1; shift ;;
    -h|--help) usage; exit 0 ;;
    --) shift; break ;;
    -*) printf 'Unknown option: %s\n' "$1" >&2; usage >&2; exit 2 ;;
    *) break ;;
  esac
done

if (($# != 2)); then
  usage >&2
  exit 2
fi

base_url=$1
file_glob=$2

[[ $jobs =~ ^[1-9][0-9]*$ ]] || { printf 'Jobs must be a positive integer.\n' >&2; exit 2; }
[[ $skip =~ ^[0-9]+$ ]] || { printf -- '--skip must be a non-negative integer.\n' >&2; exit 2; }
if [[ $num_files != -1 && ! $num_files =~ ^[0-9]+$ ]]; then
  printf -- '--numFiles must be a non-negative integer.\n' >&2
  exit 2
fi
[[ $base_url == http://* || $base_url == https://* ]] || {
  printf 'URL must begin with http:// or https://\n' >&2
  exit 2
}
[[ $base_url == */ ]] || base_url+=/

case "$tool" in
  auto)
    if command -v curl >/dev/null 2>&1; then tool=curl
    elif command -v wget >/dev/null 2>&1; then tool=wget
    else printf 'Neither curl nor wget is installed.\n' >&2; exit 1
    fi
    ;;
  curl|wget)
    command -v "$tool" >/dev/null 2>&1 || {
      printf '%s is not installed.\n' "$tool" >&2
      exit 1
    }
    ;;
  *) printf -- '--tool must be auto, curl, or wget.\n' >&2; exit 2 ;;
esac

mkdir -p "$output_dir"
tmp_dir=$(mktemp -d "${TMPDIR:-/tmp}/url-glob.XXXXXX")
trap 'rm -rf -- "$tmp_dir"' EXIT INT TERM
listing=$tmp_dir/listing.html
matches=$tmp_dir/matches.tsv
selected=$tmp_dir/selected.tsv

if [[ $tool == curl ]]; then
  curl --fail --silent --show-error --location --retry 5 --retry-all-errors \
    --output "$listing" "$base_url"
else
  wget --quiet --tries=5 --output-document="$listing" "$base_url"
fi

# Extract href values and decode percent escapes for matching/local filenames.
# The first column retains the server-provided href for the actual HTTP request.
perl -ne '
  while (/href\s*=\s*["\x27]([^"\x27]+)["\x27]/ig) {
    $href = $1;
    $href =~ s/&amp;/&/g;
    next if $href =~ m{^(?:\.\.?/|\?|#)$};
    ($name = $href) =~ s{[?#].*$}{};
    $name =~ s{.*/}{};
    $name =~ s/%([0-9A-Fa-f]{2})/chr(hex($1))/eg;
    next if $name eq "";
    print "$href\t$name\n";
  }
' "$listing" | LC_ALL=C sort -u > "$matches.all"

: > "$matches"
while IFS=$'\t' read -r href name; do
  case "$name" in
    $file_glob) printf '%s\t%s\n' "$href" "$name" >> "$matches" ;;
  esac
done < "$matches.all"

count=$(wc -l < "$matches" | tr -d ' ')
if ((count == 0)); then
  printf 'No filenames matched glob %q at %s\n' "$file_glob" "$base_url" >&2
  exit 3
fi

awk -v skip="$skip" -v limit="$num_files" '
  NR > skip && (limit < 0 || NR <= skip + limit)
' "$matches" > "$selected"

selected_count=$(wc -l < "$selected" | tr -d ' ')
printf 'Matched %s file(s); selected %s after skipping %s; downloader=%s; parallel jobs=%s\n' \
  "$count" "$selected_count" "$skip" "$tool" "$jobs"

if ((dry_run)); then
  cut -f2- "$selected"
  exit 0
fi

if ((selected_count == 0)); then
  printf 'No files selected; nothing to download.\n'
  exit 0
fi

export URL_GLOB_BASE=$base_url URL_GLOB_OUT=$output_dir URL_GLOB_TOOL=$tool

download_one() {
  local href=$1 name=$2 url final partial
  final=$URL_GLOB_OUT/$name
  partial=$final.part

  if [[ -f $final ]]; then
    printf 'SKIP  %s\n' "$name"
    return 0
  fi

  case "$href" in
    http://*|https://*) url=$href ;;
    /*)
      url=$(printf '%s' "$URL_GLOB_BASE" | sed -E 's#(https?://[^/]+).*#\1#')$href
      ;;
    *) url=$URL_GLOB_BASE$href ;;
  esac

  printf 'GET   %s\n' "$name"
  if [[ $URL_GLOB_TOOL == curl ]]; then
    curl --fail --location --retry 5 --retry-all-errors --continue-at - \
      --output "$partial" "$url"
  else
    wget --continue --tries=5 --output-document="$partial" "$url"
  fi
  mv -- "$partial" "$final"
  printf 'DONE  %s\n' "$name"
}
export -f download_one

# Launch portable batches rather than relying on GNU-only xargs options. This
# works with the older Bash and BSD utilities shipped with macOS as well.
pids=()
failures=0
wait_for_batch() {
  local pid
  for pid in "${pids[@]}"; do
    wait "$pid" || failures=$((failures + 1))
  done
  pids=()
}

while IFS=$'\t' read -r href name; do
  download_one "$href" "$name" &
  pids+=("$!")
  if ((${#pids[@]} >= jobs)); then
    wait_for_batch
  fi
done < "$selected"
wait_for_batch

if ((failures)); then
  printf '%s download(s) failed; partial files were retained for resuming.\n' \
    "$failures" >&2
  exit 1
fi

printf 'All downloads completed in %s\n' "$output_dir"

