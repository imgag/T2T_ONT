#!/usr/bin/env python3

import json
import pandas as pd
from pathlib import Path
import re
from datetime import datetime
import argparse
from typing import Optional, Dict, Any

def extract_json_from_html(file_path: Path) -> Optional[str]:
    """Extract JSON data from HTML file.
    
    Args:
        file_path: Path to HTML file
        
    Returns:
        JSON string containing report data or None if extraction fails
    """
    try:
        with open(file_path) as f:
            content = f.read()
        
        # Find the start of the JSON data
        start_match = re.search(r'const\s+reportData\s*=\s*{', content)
        if not start_match:
            print(f"Warning: No reportData found in {file_path}")
            return None
            
        # Get the position where JSON starts
        start_pos = start_match.end() - 1  # Include the opening brace
        
        # Find matching closing brace
        brace_count = 0
        pos = start_pos
        
        while pos < len(content):
            char = content[pos]
            if char == '{':
                brace_count += 1
            elif char == '}':
                brace_count -= 1
                if brace_count == 0:
                    # Found the matching closing brace
                    json_str = content[start_pos:pos + 1]
                    # Validate JSON
                    json.loads(json_str)
                    return json_str
            pos += 1
            
        print(f"Warning: Could not find end of JSON data in {file_path}")
        return None
        
    except json.JSONDecodeError as e:
        print(f"Warning: JSON parsing error in {file_path}: {str(e)}")
        return None
    except Exception as e:
        print(f"Warning: Error reading {file_path}: {str(e)}")
        return None

def safe_get(data: Dict[str, Any], *keys: str, default: Any = None) -> Any:
    """Safely get nested dictionary values.
    
    Args:
        data: Dictionary to search
        *keys: Key sequence to follow
        default: Default value if path not found
        
    Returns:
        Value at key path or default if not found
    """
    try:
        for key in keys:
            data = data[key]
        return data
    except (KeyError, TypeError, IndexError):
        return default

def find_in_setup(setup_data: list, title: str) -> Optional[str]:
    """Find value in run setup/settings data by title.
    
    Args:
        setup_data: List of setup/settings dictionaries
        title: Title to search for
        
    Returns:
        Matching value or None if not found
    """
    try:
        return next(
            (item['value'] for item in setup_data 
             if item.get('title') == title),
            None
        )
    except (TypeError, StopIteration):
        return None

def parse_report(file_path: Path) -> Optional[pd.Series]:
    """Parse report data into a pandas Series.
    
    Args:
        file_path: Path to HTML file
        
    Returns:
        Series containing parsed metadata or None if parsing fails
    """
    json_str = extract_json_from_html(file_path)
    if not json_str:
        return None
        
    try:
        data = json.loads(json_str)
    except json.JSONDecodeError as e:
        print(f"Warning: JSON parsing error in {file_path}: {str(e)}")
        return None
        
    # Extract all metadata
    metadata = {
        # Run info
        'flow_cell_id': safe_get(data, 'flow_cell_id'),
        'run_start_time': pd.to_datetime(safe_get(data, 'run_start_time')),
        'run_end_time': pd.to_datetime(safe_get(data, 'run_end_time')),
        'run_status': safe_get(data, 'run_status'),
        'run_status_context': safe_get(data, 'run_status_additional_context'),
        'run_complete': safe_get(data, 'run_complete'),
        'estimated_n50': safe_get(data, 'estimated_n50'),
        
        # Device info
        'device_type': safe_get(data, 'header', 'device_type'),
        'device_serial': safe_get(data, 'header', 'serial'),
        'experiment_name': safe_get(data, 'header', 'experiment_name'),
        'sample_id': safe_get(data, 'header', 'sample_id'),
        'position': safe_get(data, 'header', 'position'),
        'protocol_run_id': safe_get(data, 'header', 'protocol_run_id'),
        
        # Feature flags
        'alignment_enabled': safe_get(data, 'alignment_enabled'),
        'barcoding_enabled': safe_get(data, 'barcoding_enabled'),
        'basecalling_enabled': safe_get(data, 'basecalling_enabled'),
        'duplex_enabled': safe_get(data, 'duplex_enabled'),
        
        # Output stats
        'estimated_bases': safe_get(data, 'data_output', 'estimated_bases'),
        'data_produced_bytes': safe_get(data, 'data_output', 'data_produced'),
        'reads_generated': safe_get(data, 'data_output', 'reads_generated'),
        
        # Basecalling stats
        'reads_called_percent': safe_get(data, 'basecalling', 'reads_called'),
        'reads_called_pass': safe_get(data, 'basecalling', 'reads_called_pass'),
        'reads_called_fail': safe_get(data, 'basecalling', 'reads_called_fail'),
        'reads_called_skipped': safe_get(data, 'basecalling', 'reads_called_skipped'),
        'bases_called_pass': safe_get(data, 'basecalling', 'bases_called_pass'),
        
        # Run limits
        'target_bases': safe_get(data, 'run_until', 'target_estimated_bases'),
        'target_basecalled': safe_get(data, 'run_until', 'target_basecalled_bases'),
        'pores_remaining': safe_get(data, 'run_until', 'pores_remaining'),
        'run_time_limit_seconds': safe_get(data, 'run_until', 'total_run_time'),
        'elapsed_time_seconds': safe_get(data, 'run_until', 'elapsed_time_since_start_seconds'),
        
        # Run setup
        'flow_cell_type': find_in_setup(safe_get(data, 'run_setup', default=[]), 'Flow cell type'),
        'kit_type': find_in_setup(safe_get(data, 'run_setup', default=[]), 'Kit type'),
        
        # Run settings
        'basecalling_model': find_in_setup(safe_get(data, 'run_settings', default=[]), 'Basecalling'),
        'modifications': find_in_setup(safe_get(data, 'run_settings', default=[]), 'Modifications'),
        'min_qscore': find_in_setup(safe_get(data, 'run_settings', default=[]), 'Min Q score')
    }
    
    # Calculate run duration
    if metadata['run_start_time'] and metadata['run_end_time']:
        duration = metadata['run_end_time'] - metadata['run_start_time']
        metadata['run_duration_hours'] = duration.total_seconds() / 3600
    else:
        metadata['run_duration_hours'] = None
        
    return pd.Series(metadata)

def parse_reports(report_dir: Path) -> pd.DataFrame:
    """Parse multiple report files.
    
    Args:
        report_dir: Directory containing report files
        
    Returns:
        DataFrame containing parsed data from all reports
    """
    report_files = list(report_dir.glob('*.html'))
    
    if not report_files:
        raise ValueError(f"No HTML files found in {report_dir}")
        
    # Parse each file
    parsed_data = []
    for file_path in report_files:
        print(f"\nProcessing file: {file_path.name}")
        report_data = parse_report(file_path)
        if report_data is not None:
            parsed_data.append(report_data)
            
    if not parsed_data:
        raise ValueError("No valid reports could be parsed")
        
    print(f"\nSuccessfully parsed {len(parsed_data)} of {len(report_files)} files")
    
    # Combine all reports into a single DataFrame
    return pd.DataFrame(parsed_data)

def main():
    parser = argparse.ArgumentParser(description='Parse MinKNOW report HTML files and extract metadata into a CSV file.')
    parser.add_argument('report_directory', help='Directory containing MinKNOW HTML report files')
    args = parser.parse_args()
    
    report_dir = Path(args.report_directory)
    
    # Parse reports and save as CSV
    reports = parse_reports(report_dir)
    output_file = report_dir / 'run_summary.csv'
    reports.to_csv(output_file, index=False)
    print(f"Wrote {output_file}")

if __name__ == '__main__':
    main()