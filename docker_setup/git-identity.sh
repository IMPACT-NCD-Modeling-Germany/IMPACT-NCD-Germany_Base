#!/usr/bin/env bash
# s6 cont-init.d script: configures git identity + SSH for user 'rstudio'
set -euo pipefail

log() { echo "[git-init] $*" >&2; }

RSTUDIO_USER="${RSTUDIO_USER:-rstudio}"
RSTUDIO_HOME="/home/${RSTUDIO_USER}"
GIT_KNOWN_HOSTS="${GIT_KNOWN_HOSTS:-/etc/ssh/ssh_known_hosts}"
GIT_KEY_PATH="${GIT_KEY_PATH:-}"

# Helper to run a command as rstudio
as_rstudio() { su -s /bin/bash - "$RSTUDIO_USER" -c "$*"; }

log "start; USERID=${USERID:-<unset>} GROUPID=${GROUPID:-<unset>} GIT_KEY_PATH=${GIT_KEY_PATH:-<unset>}"

# 1) If key path not provided, auto-detect a single key under /keys
if [[ -z "${GIT_KEY_PATH}" && -d /keys ]]; then
  mapfile -t CANDS < <(find /keys -maxdepth 1 -type f ! -name "*.pub")
  if [[ "${#CANDS[@]}" -eq 1 ]]; then
    GIT_KEY_PATH="${CANDS[0]}"
    log "auto-detected key: ${GIT_KEY_PATH}"
  else
    log "no single key auto-detected in /keys (count=${#CANDS[@]})"
  fi
fi

# 1a) Copy SSH key to writable location and fix permissions
if [[ -n "${GIT_KEY_PATH}" && -f "${GIT_KEY_PATH}" ]]; then
  # Create a writable copy of the SSH key in rstudio home
  WRITABLE_KEY_PATH="${RSTUDIO_HOME}/.ssh/$(basename "${GIT_KEY_PATH}")"
  log "copying SSH key to writable location: ${WRITABLE_KEY_PATH}"
  
  # Ensure .ssh directory exists with proper permissions
  as_rstudio "mkdir -p ${RSTUDIO_HOME}/.ssh && chmod 700 ${RSTUDIO_HOME}/.ssh"
  
  # Copy the key and fix permissions
  cp "${GIT_KEY_PATH}" "${WRITABLE_KEY_PATH}" || { log "failed to copy SSH key"; GIT_KEY_PATH=""; }
  
  if [[ -f "${WRITABLE_KEY_PATH}" ]]; then
    chmod 600 "${WRITABLE_KEY_PATH}"
    chown "${RSTUDIO_USER}:${RSTUDIO_USER}" "${WRITABLE_KEY_PATH}" 2>/dev/null || true
    
    # Copy public key if it exists
    if [[ -f "${GIT_KEY_PATH}.pub" ]]; then
      cp "${GIT_KEY_PATH}.pub" "${WRITABLE_KEY_PATH}.pub" 2>/dev/null || true
      if [[ -f "${WRITABLE_KEY_PATH}.pub" ]]; then
        chmod 644 "${WRITABLE_KEY_PATH}.pub"
        chown "${RSTUDIO_USER}:${RSTUDIO_USER}" "${WRITABLE_KEY_PATH}.pub" 2>/dev/null || true
      fi
    fi
    
    # Update GIT_KEY_PATH to point to the writable copy
    GIT_KEY_PATH="${WRITABLE_KEY_PATH}"
    log "SSH key copied and permissions fixed: ${GIT_KEY_PATH}"
  else
    log "failed to create writable copy of SSH key"
    GIT_KEY_PATH=""
  fi
fi

# 2) Ensure GitHub host key is present (if not mounted)
if ! grep -qsE 'github\.com' "$GIT_KNOWN_HOSTS"; then
  log "adding github.com to $GIT_KNOWN_HOSTS"
  mkdir -p "$(dirname "$GIT_KNOWN_HOSTS")"
  ssh-keyscan github.com >> "$GIT_KNOWN_HOSTS" 2>/dev/null || true
fi

# 3) Make GUI (git) use the mounted key for all repos
if [[ -n "${GIT_KEY_PATH}" && -r "${GIT_KEY_PATH}" ]]; then
  as_rstudio "git config --global core.sshCommand 'ssh -i ${GIT_KEY_PATH} -o IdentitiesOnly=yes -o UserKnownHostsFile=${GIT_KNOWN_HOSTS} -o StrictHostKeyChecking=yes'"
  log "set core.sshCommand to use ${GIT_KEY_PATH}"
else
  log "no readable key at GIT_KEY_PATH; skipping core.sshCommand"
fi

# 4) Set user.name/email only if missing
NEED_NAME=1; NEED_MAIL=1
as_rstudio "git config --global --get user.name"  >/dev/null 2>&1 && NEED_NAME=0
as_rstudio "git config --global --get user.email" >/dev/null 2>&1 && NEED_MAIL=0

if [[ $NEED_NAME -eq 1 || $NEED_MAIL -eq 1 ]]; then
  GH_USER=""
  if [[ -n "${GIT_KEY_PATH}" && -r "${GIT_KEY_PATH}" ]]; then
    log "attempting to connect to GitHub with key: ${GIT_KEY_PATH}"
    
    # SSH banner prints to stderr; do not fail on nonzero exit
    # Use short timeouts suitable for container startup
    BANNER=$(ssh -T -i "${GIT_KEY_PATH}" \
      -o IdentitiesOnly=yes \
      -o UserKnownHostsFile="${GIT_KNOWN_HOSTS}" \
      -o StrictHostKeyChecking=yes \
      -o ConnectTimeout=3 \
      -o ServerAliveInterval=1 \
      -o ServerAliveCountMax=2 \
      git@github.com 2>&1 || true)
    
    log "SSH banner output: '$BANNER'"
    
    # Check if we got a connection error or empty banner
    if [[ -z "$BANNER" ]]; then
      log "no SSH banner received - network may not be available during container startup"
    elif [[ "$BANNER" =~ "Connection refused" ]] || [[ "$BANNER" =~ "No route to host" ]] || [[ "$BANNER" =~ "Network is unreachable" ]] || [[ "$BANNER" =~ "Permission denied" ]]; then
      log "SSH connection failed - will rely on fallback methods"
      BANNER=""  # Clear banner to trigger fallbacks
    fi
    
    # Try multiple patterns to extract GitHub username from banner
    # Pattern 1: Standard sed approach for "Hi username! ..."
    GH_USER=$(printf '%s' "$BANNER" | sed -n 's/^Hi \([^!]*\)!.*/\1/p' 2>/dev/null | tr -d '[:space:]' || true)
    
    # Pattern 2: AWK approach if sed fails
    if [[ -z "$GH_USER" ]]; then
      GH_USER=$(printf '%s' "$BANNER" | awk '/^Hi / { 
        match($0, /^Hi ([^!]+)!/, arr); 
        gsub(/[[:space:]]/, "", arr[1]); 
        print arr[1] 
      }' 2>/dev/null || true)
    fi
    
    # Pattern 3: Simple extraction - everything between "Hi " and "!"
    if [[ -z "$GH_USER" ]]; then
      GH_USER=$(printf '%s' "$BANNER" | sed 's/.*Hi \([^!]*\)!.*/\1/' | tr -d '[:space:]' 2>/dev/null || true)
    fi
    
    # Clean up any remaining special characters but keep valid username chars
    if [[ -n "$GH_USER" ]]; then
      GH_USER=$(printf '%s' "$GH_USER" | sed 's/[^A-Za-z0-9._-]//g' 2>/dev/null || true)
    fi
    
    [[ -n "$GH_USER" ]] && log "derived GitHub user from banner: ${GH_USER}" || log "could not derive user from banner: '$BANNER'"
  fi

  # Fallback: parse after last underscore in key filename (id_ed25519_<user>)
  if [[ -z "$GH_USER" && -n "${GIT_KEY_PATH}" ]]; then
    BASE=$(basename "$GIT_KEY_PATH")
    log "attempting to extract username from key filename: $BASE"
    
    # Try to extract username from patterns like: id_ed25519_username, id_rsa_username, etc.
    # Look for the last part after the final underscore
    if [[ "$BASE" =~ _([A-Za-z0-9._-]+)$ ]]; then
      CAND="${BASH_REMATCH[1]}"
      # Only use filename-based fallback if it looks like a valid username and isn't just the key type
      if [[ "$CAND" =~ ^[A-Za-z0-9._-]+$ && "$CAND" != "ed25519" && "$CAND" != "rsa" && "$CAND" != "dsa" && "$CAND" != "ecdsa" ]]; then
        GH_USER="$CAND"
        log "fallback user from key name: ${GH_USER}"
      else
        log "key filename pattern '$CAND' does not look like a valid GitHub username"
      fi
    else
      log "key filename '$BASE' does not contain recognizable username pattern (expected: keytype_username)"
    fi
  fi

  NAME_TO_SET="${GH_USER:-rstudio}"
  MAIL_TO_SET="${GH_USER:+${GH_USER}@users.noreply.github.com}"
  MAIL_TO_SET="${MAIL_TO_SET:-rstudio@localhost}"

  [[ $NEED_NAME -eq 1 ]] && as_rstudio "git config --global user.name  '${NAME_TO_SET}'"
  [[ $NEED_MAIL -eq 1 ]] && as_rstudio "git config --global user.email '${MAIL_TO_SET}'"
fi

# 5) Log the result
as_rstudio "echo '[git-init] user.name:  ' \$(git config --global user.name)"
as_rstudio "echo '[git-init] user.email: ' \$(git config --global user.email)"
log "done"
