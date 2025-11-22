#!/usr/bin/env bash
set -euo pipefail

# ==================================================
#  Default: run all cases
#  Options:
#    -a        => run all cases (1~6)
#    -123      => run case 1,2,3
#    -46       => run case 4,6
# ==================================================

# Parse argument (optional)
if [[ $# -eq 0 ]]; then
    CASES=(1 2 3 4 5 6)
else
    arg="$1"
    if [[ "$arg" == "-a" ]]; then
        CASES=(1 2 3 4 5 6)
    elif [[ "$arg" =~ ^-[0-9]+$ ]]; then
        arg="${arg:1}"
        CASES=($(echo "$arg" | grep -o .))
    else
        echo "用法: ./run_cases.sh [-a | -123 | -46]"
        exit 1
    fi
fi

# ==================================================
# Make the project
# ==================================================
echo "=== Running make ==="
if ! make; then
    echo "make 失敗，請修正 compile errors"
    exit 1
fi
echo "=== make 完成 ==="
echo

ROUTER="./bin/router"
EVAL="python3 utilities/pa3_evaluator.py"
INPUT_DIR="inputs"
OUTPUT_DIR="outputs"

mkdir -p "${OUTPUT_DIR}"

if [[ ! -x "${ROUTER}" ]]; then
    echo "找不到 router 執行檔：${ROUTER}"
    exit 1
fi

echo "=== PA3 Batch Runner ==="
echo "執行 cases: ${CASES[@]}"

# ==================================================
# Run all selected cases
# ==================================================
for k in "${CASES[@]}"; do
    cap="${INPUT_DIR}/case${k}.cap"
    net="${INPUT_DIR}/case${k}.net"
    route="${OUTPUT_DIR}/case${k}.route"

    echo
    echo "---------- case${k} ----------"

    if [[ ! -f "${cap}" || ! -f "${net}" ]]; then
        echo "[case${k}] 缺少輸入檔"
        exit 1
    fi

    echo "[Router] ${cap} + ${net} -> ${route}"

    tmp_time="$(mktemp)"
    /usr/bin/time -v "${ROUTER}" --cap "${cap}" --net "${net}" --out "${route}" \
        2> "${tmp_time}"

    # Parse numbers
    elapsed=$(grep -F "Elapsed (wall clock) time" "${tmp_time}" | awk -F': ' '{print $2}')
    maxrss_kb=$(grep -F "Maximum resident set size" "${tmp_time}" | awk -F': ' '{print $2}')

    if [[ -n "${maxrss_kb}" ]]; then
        maxrss_mb=$(awk "BEGIN {printf \"%.2f\", ${maxrss_kb}/1024.0}")
    else
        maxrss_mb="N/A"
    fi

    rm -f "${tmp_time}"

    echo "[case${k}] Time: ${elapsed}"
    echo "[case${k}] Max RSS: ${maxrss_mb} MB"

    echo "[Evaluator]"
    ${EVAL} "${cap}" "${net}" "${route}" || {
        echo "[case${k}] evaluator 回傳非 0"
        exit 1
    }
done

echo
echo "=== All selected cases finished ==="
