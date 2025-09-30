#!/usr/bin/env bash
# Runs at container start (s6 cont-init.d). Purpose:
# - Detect user's SSH key (mounted at /keys/…)
# - Ensure GitHub host key present
# - Derive GitHub username from SSH banner (if possible)
# - Set git user.name/email ONLY if missing
# - Make GUI (git) use that key automatically

set -euo pipefail

RSTUDIO_USER="${RSTUDIO_USER:-rstudio}"
RSTUDIO_HOME="/home/${RSTUDIO_USER}"
GIT_KNOWN_HOSTS="${GIT_KNOWN_HOSTS:-/etc/ssh/ssh_known_hosts}"
GIT_KEY_PATH="${GIT_KEY_PATH:-}"

# 1) Detect key if not provided via env
if [[ -z "${GIT_KEY_PATH}" && -d /keys ]]; then
  # Pick the only non-.pub file if exactly one exists
  mapfile -t CANDS < <(find /keys -maxdepth 1 -type f ! -name "*.pub")
  if [[ "${#CANDS[@]}" -eq 1 ]]; then
    GIT_KEY_PATH="${CANDS[0]}"
  fi
fi

# 2) Ensure GitHub host key exists (if not already mounted/populated)
if ! grep -qE 'github\.com' "$GIT_KNOWN_HOSTS" 2>/dev/null; then
  mkdir -p "$(dirname "$GIT_KNOWN_HOSTS")"
  ssh-keyscan github.com >> "$GIT_KNOWN_HOSTS" 2>/dev/null || true
fi

# Helper to run git config for the rstudio user
as_rstudio() {
  su -s /bin/bash - "$RSTUDIO_USER" -c "$*"
}

# 3) Make Git use the mounted SSH key automatically (global for rstudio user)
if [[ -n "${GIT_KEY_PATH}" && -r "${GIT_KEY_PATH}" ]]; then
  as_rstudio "git config --global core.sshCommand 'ssh -i ${GIT_KEY_PATH} -o IdentitiesOnly=yes -o UserKnownHostsFile=${GIT_KNOWN_HOSTS} -o StrictHostKeyChecking=yes -o StrictModes=no'"
fi

# 4) If user.name/email not set, try to derive from GitHub SSH banner
NEED_NAME=1; NEED_MAIL=1
as_rstudio "git config --global --get user.name" >/dev/null 2>&1 && NEED_NAME=0
as_rstudio "git config --global --get user.email" >/dev/null 2>&1 && NEED_MAIL=0

if [[ $NEED_NAME -eq 1 || $NEED_MAIL -eq 1 ]]; then
  GH_USER=""
  if [[ -n "${GIT_KEY_PATH}" && -r "${GIT_KEY_PATH}" ]]; then
    BANNER=$(ssh -T -i "${GIT_KEY_PATH}" \
      -o IdentitiesOnly=yes \
      -o UserKnownHostsFile="${GIT_KNOWN_HOSTS}" \
      -o StrictHostKeyChecking=yes \
      -o StrictModes=no \
      git@github.com 2>&1 || true)
    # Parse: "Hi <username>! You've successfully authenticated..."
    GH_USER=$(printf '%s' "$BANNER" | sed -n 's/^Hi \([^!]*\)!.*/\1/p')
  fi

  NAME_TO_SET="${GH_USER:-${RSTUDIO_USER}}"
  MAIL_TO_SET="${GH_USER:+${GH_USER}@users.noreply.github.com}"
  MAIL_TO_SET="${MAIL_TO_SET:-${RSTUDIO_USER}@localhost}"

  [[ $NEED_NAME -eq 1 ]] && as_rstudio "git config --global user.name  '${NAME_TO_SET}'"
  [[ $NEED_MAIL -eq 1 ]] && as_rstudio "git config --global user.email '${MAIL_TO_SET}'"
fi

# Done. s6 will continue to start RStudio Server.
