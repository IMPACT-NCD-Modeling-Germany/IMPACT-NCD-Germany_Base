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
        
        # Use rsync to sync files, excluding certain directories/files
        # IMPORTANT: Exclude Rpackage/ completely except R source code subfolder
        rsync -av \
            --delete \
            --exclude='docker_setup/' \
            --exclude='Rpackage/' \
            --include='Rpackage/*/R/***' \
            --exclude='*.log' \
            --exclude='.Rhistory' \
            --exclude='.RData' \
            --exclude='outputs/' \
            --exclude='inputs/synthpop/' \
            --exclude='.vscode/' \
            --exclude='*.tmp' \
            --exclude='*.cache' \
            "$REPO_SYNC_PATH/" "$CONTAINER_REPO_PATH/"
        
        # Separately sync only the R source code folders in Rpackage
        if [[ -d "$REPO_SYNC_PATH/Rpackage" ]]; then
            log "Syncing Rpackage R source folders..."
            for pkg_dir in "$REPO_SYNC_PATH/Rpackage/"*/; do
                if [[ -d "${pkg_dir}R" ]]; then
                    pkg_name=$(basename "$pkg_dir")
                    mkdir -p "$CONTAINER_REPO_PATH/Rpackage/$pkg_name"
                    rsync -av "${pkg_dir}R/" "$CONTAINER_REPO_PATH/Rpackage/$pkg_name/R/"
                fi
            done
        fi
        
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

# Configure RStudio to automatically open the R project
log "Configuring RStudio to load R project automatically..."

# Find the .Rproj file in the repository
RPROJ_FILE=$(find "$CONTAINER_REPO_PATH" -maxdepth 1 -name "*.Rproj" | head -1)

if [[ -n "$RPROJ_FILE" ]]; then
    RPROJ_NAME=$(basename "$RPROJ_FILE")
    log "Found R project file: $RPROJ_NAME"
    
    # Ensure rstudio user directories exist
    mkdir -p /home/rstudio/.rstudio/monitored/user-settings
    mkdir -p /home/rstudio/.config/rstudio
    
    # Set RStudio user preferences to open the project
    cat > /home/rstudio/.rstudio/monitored/user-settings/user-settings << EOF
{
    "initial_working_directory": "$CONTAINER_REPO_PATH",
    "default_project_location": "$CONTAINER_REPO_PATH",
    "restore_last_project": true,
    "restore_project_required": true
}
EOF

    # Create RStudio desktop file to automatically open project
    cat > /home/rstudio/.config/rstudio/rstudio-prefs.json << EOF
{
    "initial_working_directory": "$CONTAINER_REPO_PATH",
    "default_project_location": "$CONTAINER_REPO_PATH",
    "restore_last_project": true
}
EOF

    # Set the project as the last opened project
    mkdir -p /home/rstudio/.rstudio/projects_settings
    echo "$RPROJ_FILE" > /home/rstudio/.rstudio/projects_settings/last-project-path
    
    # Create project state file
    cat > /home/rstudio/.rstudio/projects_settings/project-settings << EOF
{
    "project_path": "$RPROJ_FILE"
}
EOF

    # Set proper ownership
    chown -R rstudio:rstudio /home/rstudio/.rstudio /home/rstudio/.config
    
    log "RStudio configured to open project: $RPROJ_NAME in $CONTAINER_REPO_PATH"
else
    log "No .Rproj file found in $CONTAINER_REPO_PATH - RStudio will start in repository directory"
    
    # At minimum, set the working directory
    mkdir -p /home/rstudio/.rstudio/monitored/user-settings
    cat > /home/rstudio/.rstudio/monitored/user-settings/user-settings << EOF
{
    "initial_working_directory": "$CONTAINER_REPO_PATH"
}
EOF
    chown -R rstudio:rstudio /home/rstudio/.rstudio
fi

# Setup bidirectional sync if sync is enabled and host repo is mounted
setup_bidirectional_sync() {
    if [[ "$SYNC_ENABLED" == "true" && -d "$REPO_SYNC_PATH" && "$REPO_SYNC_PATH" != "$CONTAINER_REPO_PATH" ]]; then
        log "Setting up bidirectional sync..."
        
        # Create sync back function
        sync_back_to_host() {
            log "Syncing changes back to host..."
            rsync -av \
                --delete \
                --exclude='docker_setup/' \
                --exclude='Rpackage/' \
                --include='Rpackage/*/R/***' \
                --exclude='*.log' \
                --exclude='.Rhistory' \
                --exclude='.RData' \
                --exclude='outputs/' \
                --exclude='inputs/synthpop/' \
                --exclude='.vscode/' \
                --exclude='*.tmp' \
                --exclude='*.cache' \
                --exclude='.rstudio/' \
                "$CONTAINER_REPO_PATH/" "$REPO_SYNC_PATH/"
            
            # Separately sync only the R source code folders back to host
            if [[ -d "$CONTAINER_REPO_PATH/Rpackage" ]]; then
                for pkg_dir in "$CONTAINER_REPO_PATH/Rpackage/"*/; do
                    if [[ -d "${pkg_dir}R" ]]; then
                        pkg_name=$(basename "$pkg_dir")
                        mkdir -p "$REPO_SYNC_PATH/Rpackage/$pkg_name"
                        rsync -av "${pkg_dir}R/" "$REPO_SYNC_PATH/Rpackage/$pkg_name/R/"
                    fi
                done
            fi
        }
        
        # Setup signal handler for graceful shutdown
        trap 'log "Container stopping, syncing back to host..."; sync_back_to_host; exit 0' SIGTERM SIGINT
        
        # Start background sync process
        (
            # Wait a bit for container to fully start
            sleep 30
            
            while true; do
                # Sync back every 10 seconds
                sleep 10
                
                # Only sync if there are actual changes (avoid unnecessary I/O)
                if [[ -n "$(find "$CONTAINER_REPO_PATH" -newer /tmp/last_sync_back 2>/dev/null)" ]]; then
                    sync_back_to_host
                    touch /tmp/last_sync_back
                fi
            done
        ) &
        
        # Store background process PID for cleanup
        SYNC_BACK_PID=$!
        echo $SYNC_BACK_PID > /tmp/sync_back.pid
        
        log "Bidirectional sync started (PID: $SYNC_BACK_PID)"
        log "Changes in container will sync back to host every 2 minutes"
    fi
}

# Initialize timestamp for change detection
touch /tmp/last_sync_back

# Setup bidirectional sync
setup_bidirectional_sync

log "Starting original entrypoint..."

# Execute the original RStudio entrypoint
exec /init