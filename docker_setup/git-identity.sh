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

# Debug environment information
log "start; USERID=${USERID:-<unset>} GROUPID=${GROUPID:-<unset>} GIT_KEY_PATH=${GIT_KEY_PATH:-<unset>}"
log "environment: $(uname -a)"
log "hostname: $(hostname)"
log "network interfaces: $(ip addr show 2>/dev/null | grep 'inet ' | awk '{print $2}' | tr '\n' ' ' || echo 'unavailable')"
log "mounted filesystems with /keys or ssh:"
mount | grep -E "(keys|ssh)" || log "no key-related mounts found"
log "available directories for SSH keys:"
for dir in /keys /root/.ssh /home/rstudio/.ssh /tmp/ssh /var/ssh; do
  if [[ -d "$dir" ]]; then
    log "  $dir exists ($(ls -la "$dir" 2>/dev/null | wc -l) items)"
  fi
done

# 1) If key path not provided, auto-detect a single key under /keys or common locations
if [[ -z "${GIT_KEY_PATH}" ]]; then
  log "searching for SSH keys in multiple locations..."
  
  # Search locations in order of preference
  SEARCH_PATHS=(
    "/keys"                    # Docker mount from Windows
    "/root/.ssh"              # Common container location
    "/home/rstudio/.ssh"      # User SSH directory
    "/tmp/ssh"                # Alternative mount point
    "/var/ssh"                # Another alternative
  )
  
  for search_path in "${SEARCH_PATHS[@]}"; do
    if [[ -d "$search_path" ]]; then
      log "searching for keys in: $search_path"
      log "directory contents: $(ls -la "$search_path" 2>/dev/null || echo 'cannot list')"
      
      # Look for files directly in the directory
      mapfile -t CANDS < <(find "$search_path" -maxdepth 1 -type f ! -name "*.pub" ! -name "known_hosts*" ! -name "config" 2>/dev/null || true)
      
      # Also look for SSH keys that match the id_ed25519_username pattern
      mapfile -t PATTERN_CANDS < <(find "$search_path" -maxdepth 1 -type f -name "id_ed25519_*" 2>/dev/null || true)
      CANDS+=("${PATTERN_CANDS[@]}")
      
      log "found ${#CANDS[@]} candidate files: ${CANDS[*]}"
      
      # Filter out non-key files and check if files are readable
      VALID_KEYS=()
      for candidate in "${CANDS[@]}"; do
        if [[ -f "$candidate" && -r "$candidate" ]]; then
          # Debug: show what we're checking
          log "checking candidate: $candidate"
          
          # Check if it looks like a private key by content or name
          KEY_HEADER=""
          if [[ -r "$candidate" ]]; then
            KEY_HEADER=$(head -1 "$candidate" 2>/dev/null || true)
          fi
          
          log "key header: '$KEY_HEADER'"
          
          # Check for SSH key content or naming patterns
          # Specifically look for id_ed25519_username pattern
          if [[ "$KEY_HEADER" =~ (BEGIN (RSA|DSA|EC|OPENSSH) PRIVATE KEY|BEGIN PRIVATE KEY) ]] || \
             [[ "$(basename "$candidate")" =~ ^id_ed25519_[A-Za-z0-9._-]+$ ]] || \
             [[ "$candidate" =~ /id_ed25519_[A-Za-z0-9._-]+$ ]]; then
            VALID_KEYS+=("$candidate")
            log "found valid SSH key: $candidate"
          else
            log "rejected candidate (not a recognized SSH key pattern): $candidate"
          fi
        else
          log "candidate not readable or not a file: $candidate"
        fi
      done
      
      if [[ "${#VALID_KEYS[@]}" -eq 1 ]]; then
        GIT_KEY_PATH="${VALID_KEYS[0]}"
        log "auto-detected key: ${GIT_KEY_PATH}"
        break
      elif [[ "${#VALID_KEYS[@]}" -gt 1 ]]; then
        # Multiple keys found, prefer id_ed25519_username pattern
        for key in "${VALID_KEYS[@]}"; do
          keyname=$(basename "$key")
          if [[ "$keyname" =~ ^id_ed25519_[A-Za-z0-9._-]+$ ]]; then
            GIT_KEY_PATH="$key"
            log "auto-selected preferred ed25519 key: ${GIT_KEY_PATH}"
            break
          fi
        done
        
        # If no ed25519 key, fall back to other patterns
        if [[ -z "$GIT_KEY_PATH" ]]; then
          for key in "${VALID_KEYS[@]}"; do
            keyname=$(basename "$key")
            if [[ "$keyname" =~ ^id_(rsa|dsa|ecdsa)(_.*)?$ ]]; then
              GIT_KEY_PATH="$key"
              log "auto-selected fallback key: ${GIT_KEY_PATH}"
              break
            fi
          done
        fi
        
        # If no preferred pattern, use the first one
        if [[ -z "$GIT_KEY_PATH" ]]; then
          GIT_KEY_PATH="${VALID_KEYS[0]}"
          log "multiple keys found, using first: ${GIT_KEY_PATH}"
        fi
        break
      fi
    fi
  done
  
  if [[ -z "$GIT_KEY_PATH" ]]; then
    log "no SSH keys found in any search location"
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

# 2) Ensure GitHub host key is present (handle both file and directory cases)
KNOWN_HOSTS_FILE="$GIT_KNOWN_HOSTS"

# If GIT_KNOWN_HOSTS is a directory, use a file inside it
if [[ -d "$GIT_KNOWN_HOSTS" ]]; then
  POTENTIAL_FILE="$GIT_KNOWN_HOSTS/ssh_known_hosts"
  if [[ -f "$POTENTIAL_FILE" ]]; then
    log "found existing known_hosts file: $POTENTIAL_FILE"
    KNOWN_HOSTS_FILE="$POTENTIAL_FILE"
  elif [[ -w "$GIT_KNOWN_HOSTS" ]]; then
    log "GIT_KNOWN_HOSTS is a writable directory, using file: $POTENTIAL_FILE"
    KNOWN_HOSTS_FILE="$POTENTIAL_FILE"
  else
    log "GIT_KNOWN_HOSTS directory is read-only, using persistent writable location"
    # Use a persistent location in the user's SSH directory
    mkdir -p "/home/rstudio/.ssh" 2>/dev/null || mkdir -p "/root/.ssh" 2>/dev/null || true
    if [[ -d "/home/rstudio/.ssh" ]]; then
      KNOWN_HOSTS_FILE="/home/rstudio/.ssh/known_hosts"
    else
      KNOWN_HOSTS_FILE="/root/.ssh/known_hosts"
    fi
    
    # Copy existing file if available
    if [[ -f "$POTENTIAL_FILE" ]]; then
      cp "$POTENTIAL_FILE" "$KNOWN_HOSTS_FILE" 2>/dev/null || true
      log "copied existing known_hosts to persistent location: $KNOWN_HOSTS_FILE"
    else
      log "using persistent known_hosts location: $KNOWN_HOSTS_FILE"
    fi
  fi
fi

# Create the known_hosts file if it doesn't exist
if [[ ! -f "$KNOWN_HOSTS_FILE" ]]; then
  log "creating known_hosts file: $KNOWN_HOSTS_FILE"
  mkdir -p "$(dirname "$KNOWN_HOSTS_FILE")" 2>/dev/null || true
  touch "$KNOWN_HOSTS_FILE" 2>/dev/null || {
    log "failed to create known_hosts file, using fallback location"
    # Final fallback to tmp if everything else fails
    KNOWN_HOSTS_FILE="/tmp/ssh_known_hosts"
    touch "$KNOWN_HOSTS_FILE"
    log "using temporary known_hosts file: $KNOWN_HOSTS_FILE (will not persist between container restarts)"
  }
fi

# Check if GitHub key is already present
if ! grep -qsE 'github\.com' "$KNOWN_HOSTS_FILE" 2>/dev/null; then
  log "adding github.com to $KNOWN_HOSTS_FILE"
  ssh-keyscan github.com >> "$KNOWN_HOSTS_FILE" 2>/dev/null || {
    log "failed to add github.com host key, creating basic entry"
    echo "github.com ssh-rsa AAAAB3NzaC1yc2EAAAADAQABAAABgQCj7ndNxQowgcQnjshcLrqPEiiphnt+VTTvDP6mHBL9j1aNUkY4Ue1gvwnGLVlOhGeYrnZaMgRK6+PKCUXaDbC7qtbW8gIkhL7aGCsOr/C56SJMy/BCZfxd1nWzAOxSDPgVsmerOBYfNqltV9/hWCqBywINIR+5dIg6JTJ72pcEpEjcYgXkE2YEFXV1JHnsKgbLWNlhScqb2UmyRkQyytRLtL+38TGxkxCflmO+5Z8HLUB0j8IXg7EPNT5o7pNvJgXjkLv/7sKxpwM6JWXkcEwzSY5vIJ8M6BHU+K5mZzRPEH4jm7QRVD9T/PD5VjBm8A3DfzNdOmXNO9WQGQdOUVJ8pY0pTXFDnxNWBRnYB6iTaEp2y8RzA0VQa8iJA" >> "$KNOWN_HOSTS_FILE" || true
  }
fi

# Update GIT_KNOWN_HOSTS to point to the actual file
GIT_KNOWN_HOSTS="$KNOWN_HOSTS_FILE"

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
    
    # Test SSH connectivity with progressively more permissive settings
    BANNER=""
    
    # Try 1: Strict settings (for local/direct connections)
    log "trying strict SSH connection..."
    BANNER=$(ssh -T -i "${GIT_KEY_PATH}" \
      -o IdentitiesOnly=yes \
      -o UserKnownHostsFile="${GIT_KNOWN_HOSTS}" \
      -o StrictHostKeyChecking=yes \
      -o ConnectTimeout=5 \
      -o ServerAliveInterval=2 \
      -o ServerAliveCountMax=3 \
      git@github.com 2>&1 || true)
    
    # Try 2: Relaxed host checking (for remote/proxied connections)
    if [[ -z "$BANNER" ]] || [[ "$BANNER" =~ (Connection refused|No route to host|Network is unreachable|timeout) ]]; then
      log "strict connection failed, trying with relaxed host checking..."
      BANNER=$(ssh -T -i "${GIT_KEY_PATH}" \
        -o IdentitiesOnly=yes \
        -o StrictHostKeyChecking=no \
        -o UserKnownHostsFile=/dev/null \
        -o ConnectTimeout=10 \
        -o ServerAliveInterval=3 \
        -o ServerAliveCountMax=5 \
        git@github.com 2>&1 || true)
    fi
    
    # Try 3: Very permissive (for complex network setups)
    if [[ -z "$BANNER" ]] || [[ "$BANNER" =~ (Connection refused|No route to host|Network is unreachable|timeout) ]]; then
      log "relaxed connection failed, trying very permissive settings..."
      BANNER=$(timeout 15 ssh -T -i "${GIT_KEY_PATH}" \
        -o IdentitiesOnly=yes \
        -o StrictHostKeyChecking=no \
        -o UserKnownHostsFile=/dev/null \
        -o ConnectTimeout=15 \
        -o BatchMode=yes \
        -o LogLevel=ERROR \
        git@github.com 2>&1 || true)
    fi
    
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
