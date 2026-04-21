#!/usr/bin/env bash

set -euo pipefail

SESSION_NAME="${1:-gcmtl-paper}"
REPO_DIR="${REPO_DIR:-/Users/diren/Documents/GC_MTL_repo_https}"
PAPER_DIR="${PAPER_DIR:-/Users/diren/Documents/GC_MTL}"
PROMPT_FILE="${PROMPT_FILE:-$REPO_DIR/scripts/prompts/gcmtl_paper_supervisor_prompt.txt}"
MODEL="${MODEL:-gpt-5.4}"
LOG_DIR="${LOG_DIR:-$REPO_DIR/logs/tmux}"

mkdir -p "$LOG_DIR"

if ! command -v tmux >/dev/null 2>&1; then
  echo "tmux is not installed." >&2
  exit 1
fi

if ! command -v codex >/dev/null 2>&1; then
  echo "codex CLI is not installed." >&2
  exit 1
fi

if [[ ! -d "$REPO_DIR" ]]; then
  echo "Repository directory not found: $REPO_DIR" >&2
  exit 1
fi

if [[ ! -f "$PROMPT_FILE" ]]; then
  echo "Prompt file not found: $PROMPT_FILE" >&2
  exit 1
fi

if tmux has-session -t "$SESSION_NAME" 2>/dev/null; then
  echo "tmux session already exists: $SESSION_NAME" >&2
  echo "Attach with: tmux attach -t $SESSION_NAME" >&2
  exit 1
fi

STAMP="$(date +%Y%m%d_%H%M%S)"
LOG_FILE="$LOG_DIR/${SESSION_NAME}_${STAMP}.log"

TMUX_CMD=$(cat <<EOF
cd "$REPO_DIR" && \
codex exec \
  --dangerously-bypass-approvals-and-sandbox \
  -m "$MODEL" \
  -C "$REPO_DIR" \
  --add-dir "$PAPER_DIR" \
  - < "$PROMPT_FILE" 2>&1 | tee "$LOG_FILE"
EOF
)

tmux new-session -d -s "$SESSION_NAME" "$TMUX_CMD"

echo "Started tmux session: $SESSION_NAME"
echo "Attach: tmux attach -t $SESSION_NAME"
echo "Log: $LOG_FILE"
