#!/usr/bin/env python3
"""
Script to automatically queue backup jobs for all folders in a given directory.
Tracks progress using marker files to avoid duplicate backup jobs.
"""

import os
import sys
import argparse
import subprocess
import logging
from pathlib import Path
from datetime import datetime
import json

class BackupManager:
    def __init__(self, base_folder, mode="run", email=None, dry_run=False, force=False):
        self.base_folder = Path(base_folder).resolve()
        self.mode = mode
        self.email = email
        self.dry_run = dry_run
        self.force = force
        self.marker_file = ".backup_queued"
        self.log_file = self.base_folder / "backup_queue.log"
        
        # Setup logging
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(self.log_file),
                logging.StreamHandler()
            ]
        )
        self.logger = logging.getLogger(__name__)
    
    def has_backup_marker(self, folder_path):
        """Check if backup marker file exists."""
        marker_path = folder_path / self.marker_file
        return marker_path.exists()
    
    def get_backup_status(self, folder_path):
        """Get detailed backup status from marker file."""
        marker_path = folder_path / self.marker_file
        if not marker_path.exists():
            return None
        
        try:
            with open(marker_path, 'r') as f:
                return json.load(f)
        except (json.JSONDecodeError, IOError):
            # If marker file is corrupted, treat as not backed up
            return None
    
    def create_backup_marker(self, folder_path, sge_job_id=None):
        """Create backup marker file with metadata."""
        marker_path = folder_path / self.marker_file
        marker_data = {
            "queued_timestamp": datetime.now().isoformat(),
            "mode": self.mode,
            "sge_job_id": sge_job_id,
            "status": "queued"
        }
        
        with open(marker_path, 'w') as f:
            json.dump(marker_data, f, indent=2)
    
    def queue_backup(self, folder_path):
        """Queue backup for a single folder."""
        folder_path = Path(folder_path).resolve()
        
        # Construct the backup command
        cmd = [
            "sudo", "-u", "archive-gs", "php", 
            f"/mnt/storage2/megSAP/pipeline/src/IMGAG/backup_queue.php",
            "-in", str(folder_path),
            "-mode", self.mode,
            "-email", "ahgrosc1"
        ]
        
        if self.email:
            cmd.extend(["-email", self.email])
        
        self.logger.info(f"Queuing backup for: {folder_path}")
        self.logger.info(f"Command: {' '.join(cmd)}")
        
        if self.dry_run:
            self.logger.info("DRY RUN - Would execute above command")
            return "dry_run"
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            
            # Parse SGE job ID from output
            sge_job_id = None
            for line in result.stdout.split('\n'):
                if 'SGE job id:' in line:
                    sge_job_id = line.split(':')[1].strip()
                    break
            
            self.logger.info(f"Successfully queued backup. Job ID: {sge_job_id}")
            return sge_job_id
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Failed to queue backup for {folder_path}: {e}")
            self.logger.error(f"stdout: {e.stdout}")
            self.logger.error(f"stderr: {e.stderr}")
            return None
    
    def get_all_folders(self):
        """Get all folders in the base directory."""
        folders = []
        
        if not self.base_folder.exists():
            self.logger.error(f"Base folder does not exist: {self.base_folder}")
            return folders
        
        for item in self.base_folder.iterdir():
            if item.is_dir():
                folders.append(item)
        
        return sorted(folders)
    
    def process_all_folders(self):
        """Process all folders in the base directory."""
        folders = self.get_all_folders()
        
        if not folders:
            self.logger.info("No folders found.")
            return
        
        self.logger.info(f"Found {len(folders)} folders")
        
        queued_count = 0
        skipped_count = 0
        failed_count = 0
        
        for folder in folders:
            self.logger.info(f"Processing folder: {folder.name}")
            
            # Check if backup already queued/completed
            if not self.force and self.has_backup_marker(folder):
                backup_status = self.get_backup_status(folder)
                if backup_status:
                    self.logger.info(f"Skipping {folder.name} - already processed on {backup_status.get('queued_timestamp', 'unknown date')}")
                else:
                    self.logger.info(f"Skipping {folder.name} - backup marker exists")
                skipped_count += 1
                continue
            
            # Queue the backup
            sge_job_id = self.queue_backup(folder)
            
            if sge_job_id is not None:
                # Only create marker file if not in dry run mode
                if not self.dry_run:
                    self.create_backup_marker(folder, sge_job_id)
                queued_count += 1
            else:
                failed_count += 1
        
        # Summary
        self.logger.info(f"Processing complete:")
        self.logger.info(f"  Queued: {queued_count}")
        self.logger.info(f"  Skipped: {skipped_count}")
        self.logger.info(f"  Failed: {failed_count}")
    
    def list_status(self):
        """List the backup status of all folders."""
        folders = self.get_all_folders()
        
        print(f"{'Folder':<40} {'Status':<15} {'Queued Date':<20} {'Job ID':<10}")
        print("-" * 90)
        
        for folder in folders:
            status_info = self.get_backup_status(folder)
            if status_info:
                status = status_info.get('status', 'queued')
                queued_date = status_info.get('queued_timestamp', 'unknown')[:19]
                job_id = status_info.get('sge_job_id', 'N/A')
            else:
                status = 'not_queued'
                queued_date = 'N/A'
                job_id = 'N/A'
            
            print(f"{folder.name:<40} {status:<15} {queued_date:<20} {job_id:<10}")


def main():
    parser = argparse.ArgumentParser(
        description="Queue backup jobs for all folders in a directory"
    )
    parser.add_argument(
        "base_folder",
        help="Base directory containing folders"
    )
    parser.add_argument(
        "--mode",
        choices=["run", "project", "user"],
        default="run",
        help="Backup mode (default: run)"
    )
    parser.add_argument(
        "--email",
        help="Email for SGE job notifications"
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be done without actually queuing jobs"
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Force backup even if marker file exists"
    )
    parser.add_argument(
        "--list-status",
        action="store_true",
        help="List backup status of all folders and exit"
    )
    
    args = parser.parse_args()
    
    # Validate base folder
    if not os.path.exists(args.base_folder):
        print(f"Error: Base folder '{args.base_folder}' does not exist")
        sys.exit(1)
    
    # Create backup manager
    backup_manager = BackupManager(
        base_folder=args.base_folder,
        mode=args.mode,
        email=args.email,
        dry_run=args.dry_run,
        force=args.force
    )
    
    if args.list_status:
        backup_manager.list_status()
    else:
        backup_manager.process_all_folders()


if __name__ == "__main__":
    main()