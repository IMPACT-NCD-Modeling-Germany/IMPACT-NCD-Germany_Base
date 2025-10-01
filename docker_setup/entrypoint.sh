#!/bin/bash
# -----------------------------------------------------------------------------
# entrypoint.sh - Docker Container Repository Sync and Entrypoint
# -----------------------------------------------------------------------------
#
# PURPOSE:
# This entrypoint script handles syncing of repository changes from the host
# to the container at startup. Since the Docker image contains a snapshot of
# the repository at build time, this script ensures the container has the
# latest changes when it starts.
#
# HOW IT WORKS:
# 1. Checks for a mounted repository at REPO_SYNC_PATH
# 2. Syncs changed files to the container's working directory
# 3. Preserves important container-specific configurations
# 4. Starts the original container entrypoint
#
# ENVIRONMENT VARIABLES:
# - REPO_SYNC_PATH: Path where the host repository is mounted (default: /host-repo)
# - SYNC_ENABLED: Enable/disable syncing (default: true)
#
# -----------------------------------------------------------------------------

# Configuration
REPO_SYNC_PATH="${REPO_SYNC_PATH:-/host-repo}"
SYNC_ENABLED="${SYNC_ENABLED:-true}"

# Get repository name from environment (set during build) or use default
source /etc/environment 2>/dev/null || true
CONTAINER_REPO_NAME="${CONTAINER_REPO_NAME:-IMPACT-NCD-Germany_Base}"
CONTAINER_REPO_PATH="/home/rstudio/$CONTAINER_REPO_NAME"

log() {
    echo "[ENTRYPOINT] $(date '+%Y-%m-%d %H:%M:%S') - $*"
}

sync_repository() {
    log "Starting repository sync..."
    
    # Check if mounted repository exists and is a git repository
    if [[ ! -d "$REPO_SYNC_PATH" ]]; then
        log "No repository mounted at $REPO_SYNC_PATH - using image snapshot"
        return 0
    fi
    
    if [[ ! -d "$REPO_SYNC_PATH/.git" ]]; then
        log "Mounted path $REPO_SYNC_PATH is not a git repository - using image snapshot"
        return 0
    fi
    
    log "Repository mounted at $REPO_SYNC_PATH - using mounted version directly"
    
    # Option 1: Direct mount (no sync needed)
    if [[ "$REPO_SYNC_PATH" != "$CONTAINER_REPO_PATH" ]]; then
        log "Syncing from $REPO_SYNC_PATH to $CONTAINER_REPO_PATH"
        
        # Ensure container repository directory exists
        #mkdir -p "$CONTAINER_REPO_PATH"
        
        # Use rsync to sync files, excluding certain directories/files
        rsync -av \
            --delete \
            --exclude='docker_setup/' \
            --exclude='*.log' \
            --exclude='.Rhistory' \
            --exclude='.RData' \
            --exclude='outputs/' \
            --exclude='inputs/synthpop/' \
            --exclude='.vscode/' \
            --exclude='*.tmp' \
            --exclude='*.cache' \
            "$REPO_SYNC_PATH/" "$CONTAINER_REPO_PATH/"
        
        # Ensure correct ownership
        chown -R rstudio:rstudio "$CONTAINER_REPO_PATH"
        chmod -R u+w "$CONTAINER_REPO_PATH"
        
        log "Repository sync completed successfully"
    else
        log "Repository already mounted at correct location"
    fi
}

# Main execution
log "Container starting..."

# Perform repository sync if enabled
if [[ "$SYNC_ENABLED" == "true" ]]; then
    sync_repository
else
    log "Repository sync disabled (SYNC_ENABLED=$SYNC_ENABLED)"
fi

# Change to the repository directory
cd "$CONTAINER_REPO_PATH" || {
    log "Warning: Could not change to $CONTAINER_REPO_PATH"
}

log "Starting original entrypoint..."

# Execute the original RStudio entrypoint
exec /init