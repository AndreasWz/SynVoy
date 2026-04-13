#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RUNTIME_DIR="$ROOT_DIR/local_data/ollama_runtime"
OLLAMA_BIN="$RUNTIME_DIR/bin/ollama"
OLLAMA_HOME_DIR="$ROOT_DIR/local_data/ollama_home"
OLLAMA_MODELS_DIR="$ROOT_DIR/local_data/ollama_models"
OLLAMA_LOG="$ROOT_DIR/local_data/ollama_server.log"
OLLAMA_HOST_URL="${OLLAMA_HOST:-http://127.0.0.1:11434}"

usage() {
    cat <<EOF
Usage: $(basename "$0") <start|status|pull|list> [args]

Commands:
  start         Start the local Ollama server in the background
  status        Check whether the local Ollama API is reachable
  pull [model]  Pull a model tag (default: gemma4:e4b)
  list          List locally available models
EOF
}

require_runtime() {
    if [[ ! -x "$OLLAMA_BIN" ]]; then
        echo "Missing local Ollama runtime at $OLLAMA_BIN" >&2
        echo "Download it into local_data/ollama_runtime first." >&2
        exit 1
    fi
}

ollama_cmd() {
    HOME="$OLLAMA_HOME_DIR" \
    OLLAMA_MODELS="$OLLAMA_MODELS_DIR" \
    "$OLLAMA_BIN" "$@"
}

ensure_dirs() {
    mkdir -p "$OLLAMA_HOME_DIR" "$OLLAMA_MODELS_DIR"
}

is_running() {
    curl -sf "$OLLAMA_HOST_URL/api/version" >/dev/null
}

start_server() {
    require_runtime
    ensure_dirs

    if is_running; then
        echo "Ollama is already running at $OLLAMA_HOST_URL"
        return 0
    fi

    nohup env \
        HOME="$OLLAMA_HOME_DIR" \
        OLLAMA_MODELS="$OLLAMA_MODELS_DIR" \
        "$OLLAMA_BIN" serve \
        >"$OLLAMA_LOG" 2>&1 &
    disown  # detach from shell so VSCode/terminal crashes don't kill it

    sleep 2
    if is_running; then
        echo "Started Ollama at $OLLAMA_HOST_URL"
        echo "Log: $OLLAMA_LOG"
        return 0
    fi

    echo "Ollama did not start successfully. Check $OLLAMA_LOG" >&2
    return 1
}

show_status() {
    if is_running; then
        echo "Ollama is running at $OLLAMA_HOST_URL"
        curl -sf "$OLLAMA_HOST_URL/api/version"
    else
        echo "Ollama is not running at $OLLAMA_HOST_URL" >&2
        return 1
    fi
}

pull_model() {
    require_runtime
    ensure_dirs
    local model="${1:-gemma4:e4b}"

    if ! is_running; then
        echo "Ollama is not running. Starting it first..." >&2
        start_server
    fi

    ollama_cmd pull "$model"
}

list_models() {
    require_runtime
    ensure_dirs
    ollama_cmd list
}

main() {
    local cmd="${1:-}"
    case "$cmd" in
        start)
            start_server
            ;;
        status)
            show_status
            ;;
        pull)
            shift || true
            pull_model "${1:-gemma4:e4b}"
            ;;
        list)
            list_models
            ;;
        *)
            usage
            exit 1
            ;;
    esac
}

main "$@"
