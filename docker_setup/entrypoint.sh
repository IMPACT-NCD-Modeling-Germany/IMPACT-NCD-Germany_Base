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
        # Preserve compiled R package files but sync source code
        rsync -av \
            --delete \
            --exclude='docker_setup/' \
            --exclude='Rpackage/*/DESCRIPTION' \
            --exclude='Rpackage/*/NAMESPACE' \
            --exclude='Rpackage/*/man/' \
            --exclude='Rpackage/*/*.tar.gz' \
            --exclude='Rpackage/*/src/*.o' \
            --exclude='Rpackage/*/src/*.so' \
            --exclude='Rpackage/*/src/*.dll' \
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
            log "DEBUG: Container path: $CONTAINER_REPO_PATH"
            log "DEBUG: Host path: $REPO_SYNC_PATH"
            log "DEBUG: Checking if paths exist..."
            ls -la "$CONTAINER_REPO_PATH" >/dev/null 2>&1 && log "DEBUG: Container path exists" || log "DEBUG: Container path MISSING"
            ls -la "$REPO_SYNC_PATH" >/dev/null 2>&1 && log "DEBUG: Host path exists" || log "DEBUG: Host path MISSING"
            
            # Show what we're about to sync
            log "DEBUG: Files in container analysis folder:"
            ls -la "$CONTAINER_REPO_PATH/analysis/" 2>/dev/null | head -5 || log "DEBUG: No analysis folder found"
            
            # Main sync back - same exclusions as forward sync but preserve compiled R package
            log "DEBUG: Running incremental rsync from $CONTAINER_REPO_PATH/ to $REPO_SYNC_PATH/"
            rsync -av \
                --update \
                --inplace \
                --exclude='docker_setup/' \
                --exclude='Rpackage/*/DESCRIPTION' \
                --exclude='Rpackage/*/NAMESPACE' \
                --exclude='Rpackage/*/man/' \
                --exclude='Rpackage/*/*.tar.gz' \
                --exclude='Rpackage/*/src/*.o' \
                --exclude='Rpackage/*/src/*.so' \
                --exclude='Rpackage/*/src/*.dll' \
                --exclude='.Rproj.user/' \
                --exclude='.git/index' \
                --exclude='.git/logs/' \
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
            
            RSYNC_EXIT_CODE=$?
            log "DEBUG: rsync exit code: $RSYNC_EXIT_CODE"
            
            if [ $RSYNC_EXIT_CODE -ne 0 ]; then
                log "ERROR: Main rsync back failed with exit code $RSYNC_EXIT_CODE"
            else
                log "DEBUG: rsync completed successfully"
            fi
            
            log "DEBUG: Sync back completed"
        }
        
        # Setup signal handler for graceful shutdown
        trap 'log "Container stopping, syncing back to host..."; sync_back_to_host; log "Final sync completed"; exit 0' SIGTERM SIGINT
        
        # Start background sync process
        (
            # Wait a bit for container to fully start
            sleep 30
            log "Background sync process started"
            
            while true; do
                # Sync back every 10 seconds
                sleep 3
                
                # Only sync if there are actual changes (avoid unnecessary I/O)
                # Exclude RStudio temp files and Git index from change detection
                ALL_CHANGED=$(find "$CONTAINER_REPO_PATH" -newer /tmp/last_sync_back 2>/dev/null)
                CHANGED_FILES=$(echo "$ALL_CHANGED" | \
                    grep -v '\.Rproj\.user' | \
                    grep -v '\.rstudio' | \
                    grep -v '\.git/index' | \
                    grep -v '\.git/logs')
                
                FILTERED_COUNT=$(echo "$ALL_CHANGED" | wc -l)
                FINAL_COUNT=$(echo "$CHANGED_FILES" | grep -v '^$' | wc -l)
                log "DEBUG: Found $FILTERED_COUNT total changes, filtered to $FINAL_COUNT relevant changes"
                if [[ -n "$CHANGED_FILES" ]]; then
                    log "Changes detected, syncing back to host..."
                    log "DEBUG: Number of changed files: $(echo "$CHANGED_FILES" | grep -v '^$' | wc -l)"
                    echo "$CHANGED_FILES" | grep -v '^$' | head -5 | while read file; do
                        log "DEBUG: Changed: $file"
                    done
                    if [[ $(echo "$CHANGED_FILES" | grep -v '^$' | wc -l) -gt 5 ]]; then
                        log "DEBUG: ... and $(($(echo "$CHANGED_FILES" | grep -v '^$' | wc -l) - 5)) more files"
                    fi
                    
                    # Create list of changed files for rsync (convert to relative paths)
                    echo "$CHANGED_FILES" | grep -v '^$' | while read file; do
                        # Convert absolute path to relative path
                        if [[ -n "$file" && "$file" == "$CONTAINER_REPO_PATH"/* ]]; then
                            REL_PATH="${file#$CONTAINER_REPO_PATH/}"
                            # Only include if it's not just the root directory
                            if [[ -n "$REL_PATH" && "$REL_PATH" != "$CONTAINER_REPO_PATH" ]]; then
                                echo "$REL_PATH"
                            fi
                        fi
                    done > /tmp/changed_files_list
                    
                    # Only sync the changed files
                    if [[ -s /tmp/changed_files_list ]]; then
                        rsync -av \
                            --files-from=/tmp/changed_files_list \
                            --exclude='docker_setup/' \
                            --exclude='Rpackage/*/DESCRIPTION' \
                            --exclude='Rpackage/*/NAMESPACE' \
                            --exclude='Rpackage/*/man/' \
                            --exclude='Rpackage/*/*.tar.gz' \
                            --exclude='Rpackage/*/src/*.o' \
                            --exclude='Rpackage/*/src/*.so' \
                            --exclude='Rpackage/*/src/*.dll' \
                            --exclude='.Rproj.user/' \
                            --exclude='.git/index' \
                            --exclude='.git/logs/' \
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
                        
                        log "DEBUG: Synced $(cat /tmp/changed_files_list | wc -l) changed files"
                    fi
                    
                    touch /tmp/last_sync_back
                else
                    log "DEBUG: No changes detected since last sync ($(date))"
                fi
            done
        ) &
        
        # Store background process PID for cleanup
        SYNC_BACK_PID=$!
        echo $SYNC_BACK_PID > /tmp/sync_back.pid
        
        log "Bidirectional sync started (PID: $SYNC_BACK_PID)"
        log "Changes in container will sync back to host every 3 seconds"
    fi
}

# Initialize timestamp for change detection
touch /tmp/last_sync_back

# Setup bidirectional sync
setup_bidirectional_sync

log "Starting original entrypoint..."

# Execute the original RStudio entrypoint
exec /init