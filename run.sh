#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 ]]; then
    echo "用法: ./run.sh [-a | -123 | -46] <name>"
    exit 1
fi

arg="$1"
name="$2"

if [[ "$arg" == "-a" ]]; then
    CASES=(1 2 3 4 5 6)
elif [[ "$arg" =~ ^-[0-9]+$ ]]; then
    arg="${arg:1}"
    CASES=($(echo "$arg" | grep -o .))
else
    echo "用法: ./run.sh [-a | -123 | -46] <name>"
    exit 1
fi

echo "=== Running make ==="
make clean
if ! make; then
    echo "make 失敗，請修正 compile errors"
    exit 1
fi
echo

ROUTER="./bin/router"
EVAL="python3 utilities/pa3_evaluator.py"
INPUT_DIR="inputs"
OUTPUT_DIR="outputs"
EVAL_DIR="evaluate/${name}"
SUMMARY="${EVAL_DIR}/summary.txt"

mkdir -p "${OUTPUT_DIR}" "${EVAL_DIR}"

if [[ ! -x "${ROUTER}" ]]; then
    echo "找不到 router 執行檔：${ROUTER}"
    exit 1
fi

echo "=== PA3 Batch Runner ==="
echo "執行 cases: ${CASES[@]}"
echo "summary 輸出在：${SUMMARY}"

# 清空 summary
: > "${SUMMARY}"

for k in "${CASES[@]}"; do
    cap="${INPUT_DIR}/case${k}.cap"
    net="${INPUT_DIR}/case${k}.net"
    route="${OUTPUT_DIR}/case${k}.route"
    eval_out="${EVAL_DIR}/case${k}.txt"

    echo
    echo "---------- case${k} ----------"

    if [[ ! -f "${cap}" || ! -f "${net}" ]]; then
        echo "[case${k}] 缺少輸入檔"
        exit 1
    fi

    echo "[Router] ${cap} + ${net} -> ${route}"
    time "${ROUTER}" --cap "${cap}" --net "${net}" --out "${route}"

    echo "[Evaluator] → ${eval_out}"
    result="$(${EVAL} "${cap}" "${net}" "${route}")"

    # 個別輸出
    # echo "${result}" > "${eval_out}"

    # append 到 summary
    {
        echo "===== case${k} ====="
        echo "${result}"
        echo
    } >> "${SUMMARY}"
done

echo
echo "=== All selected cases finished ==="
echo "summary 已產生：${SUMMARY}"
