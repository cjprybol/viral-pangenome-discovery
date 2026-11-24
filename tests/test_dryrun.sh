#!/usr/bin/env bash
# Smoke test: every Snakefile parses, and each parameterized module resolves
# its job DAG against the example configs in config/examples/ and the tiny
# synthetic inputs in tests/data/. Nothing is executed (dry-run only), so no
# bioinformatics tools, databases, or network access are needed.
#
# `snakemake --lint` is also run for each module and its findings are
# summarized. Lint findings are reported but do NOT fail the test: the 2022
# code predates several current Snakemake style rules (log directives,
# os.path.join instead of '+', etc.). Pass --strict-lint to make them fatal.
#
# Usage: tests/test_dryrun.sh [--strict-lint]
#   SNAKEMAKE=/path/to/snakemake tests/test_dryrun.sh   # override the binary

set -euo pipefail

usage() {
    sed -n '2,13p' "$0" | sed 's/^# \{0,1\}//'
}

strict_lint=0
for arg in "$@"; do
    case "${arg}" in
        --strict-lint) strict_lint=1 ;;
        --help) usage; exit 0 ;;
        *) echo "unknown argument: ${arg}" >&2; usage >&2; exit 2 ;;
    esac
done

repo_root="$(cd "$(dirname "$0")/.." && pwd)"
cd "${repo_root}"
snakemake_bin="${SNAKEMAKE:-snakemake}"
command -v "${snakemake_bin}" >/dev/null || { echo "snakemake not found: ${snakemake_bin}" >&2; exit 2; }
echo "Using $("${snakemake_bin}" --version 2>/dev/null | head -n 1) from $(command -v "${snakemake_bin}")"

log="$(mktemp)"
trap 'rm -f "${log}"' EXIT
failures=0
lint_findings=0

# module:target(s). Target list is space-separated; empty = default rule.
modules="
sra-download:
fastq-qc:all
fastq-subset:
assemble:megahit_assemble spades_assemble
classify-reads:
map-reads:bwa_mem_align bowtie2_index_reference_fasta
ncbi-datasets:merge_fastas
download-databases:download_kraken_pluspfp8_db download_kraken_pluspfp16_db
"

while IFS= read -r entry; do
    [[ -z "${entry}" ]] && continue
    module="${entry%%:*}"
    targets="${entry#*:}"
    snakefile="snakemake/${module}.snakemake.py"
    configfile="config/examples/${module}.yaml"

    # shellcheck disable=SC2086  # targets is intentionally word-split
    if "${snakemake_bin}" --snakefile "${snakefile}" --configfile "${configfile}" \
        --dry-run --cores 1 ${targets} >"${log}" 2>&1; then
        jobs="$(awk '$1 == "total" {print $2; exit}' "${log}")"
        echo "PASS dry-run  ${module} (${jobs:-?} jobs)"
    else
        echo "FAIL dry-run  ${module}"
        sed 's/^/    /' "${log}"
        failures=$((failures + 1))
    fi

    if "${snakemake_bin}" --snakefile "${snakefile}" --configfile "${configfile}" \
        --lint >"${log}" 2>&1; then
        echo "     lint     ${module}: clean"
    else
        n="$(grep -c '^    \*' "${log}" || true)"
        if [[ "${n}" -eq 0 ]]; then
            echo "FAIL lint     ${module}: lint did not run"
            sed 's/^/    /' "${log}"
            failures=$((failures + 1))
        else
            echo "     lint     ${module}: ${n} finding(s)"
            lint_findings=$((lint_findings + n))
        fi
    fi
done <<EOF
${modules}
EOF

# Exploratory modules with hard-coded paths: check they parse only.
for module in classify-contigs map-contigs; do
    if "${snakemake_bin}" --snakefile "snakemake/${module}.snakemake.py" \
        --config reference_fasta=unused query_fasta=unused \
        --list-rules >"${log}" 2>&1; then
        echo "PASS parse    ${module} ($(grep -c . "${log}") rules)"
    else
        echo "FAIL parse    ${module}"
        sed 's/^/    /' "${log}"
        failures=$((failures + 1))
    fi
done

echo
echo "lint findings (informational): ${lint_findings}"
if [[ "${strict_lint}" -eq 1 && "${lint_findings}" -gt 0 ]]; then
    echo "--strict-lint set: treating lint findings as failures"
    failures=$((failures + 1))
fi
if [[ "${failures}" -gt 0 ]]; then
    echo "${failures} failure(s)"
    exit 1
fi
echo "all checks passed"
