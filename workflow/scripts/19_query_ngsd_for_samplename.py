import os
import pandas as pd
import mysql.connector  # Fix the import
from typing import Dict, List
import sys

# File path configuration
TSV_FILE = "doc/tables/flowcell_biological_sample.tsv"


def read_existing_data(file_path: str) -> pd.DataFrame:
    """Read existing TSV file if it exists."""
    if os.path.exists(file_path):
        return pd.read_csv(file_path, sep="\t")
    return pd.DataFrame(
        columns=["name_ngsd", "name_external", "project_name", "run_flowcell_id"]
    )


def get_database_connection():
    """Create database connection using environment variables."""
    try:
        connection = mysql.connector.connect(
            host=os.getenv("MYSQL_HOST", "localhost"),
            user=os.getenv("MYSQL_USER"),
            password=os.getenv("MYSQL_PASSWORD"),
            database=os.getenv("MYSQL_DATABASE"),
            auth_plugin="mysql_native_password",  # Add authentication method
        )
        return connection
    except mysql.connector.Error as err:
        print(f"Error connecting to database: {err}")
        sys.exit(1)


def query_database(connection) -> pd.DataFrame:
    """Query the database for flowcell and sample information."""
    query = """
    SELECT
        CONCAT(sample.name, '_', LPAD(processed_sample.process_id, 2, '0')) as name_ngsd,
        sample.name_external,
        project.name as project_name,
        sequencing_run.fcid as run_flowcell_id,
        processing_system.name_short as processing_system_name
    FROM
        processed_sample
        JOIN sample ON processed_sample.sample_id = sample.id
        JOIN project ON processed_sample.project_id = project.id
        JOIN processing_system ON processed_sample.processing_system_id = processing_system.id
        JOIN sequencing_run ON processed_sample.sequencing_run_id = sequencing_run.id
    WHERE
        project.name = "25006_1422_BEGIN_T2T_GoE"
    ORDER BY
        sample.name
    """

    df = pd.read_sql(query, connection)
    return df


def merge_data(existing_df: pd.DataFrame, new_df: pd.DataFrame) -> pd.DataFrame:
    """Merge existing and new data, preserving all columns and manual changes."""
    # Create a composite key for comparison
    existing_df["composite_key"] = existing_df.apply(
        lambda x: f"{x['name_ngsd']}_{x['run_flowcell_id']}", axis=1
    )
    new_df["composite_key"] = new_df.apply(
        lambda x: f"{x['name_ngsd']}_{x['run_flowcell_id']}", axis=1
    )

    # Keep existing records as is
    result_df = existing_df.copy()

    # Find new records from database
    new_records = new_df[~new_df["composite_key"].isin(existing_df["composite_key"])]

    # Append only new records
    if len(new_records) > 0:
        # Make sure new_records has all columns from existing_df
        for col in existing_df.columns:
            if col not in new_records.columns and col != "composite_key":
                new_records[col] = None

        result_df = pd.concat([result_df, new_records[result_df.columns]])

    # Drop the composite key column
    result_df = result_df.drop("composite_key", axis=1)

    # Sort by name_ngsd
    return result_df.sort_values("name_ngsd").reset_index(drop=True)


def main():
    # Check for required environment variables
    required_env_vars = ["MYSQL_USER", "MYSQL_PASSWORD", "MYSQL_DATABASE"]
    missing_vars = [var for var in required_env_vars if not os.getenv(var)]

    if missing_vars:
        print(f"Missing required environment variables: {', '.join(missing_vars)}")
        sys.exit(1)

    # Read existing data
    existing_df = read_existing_data(TSV_FILE)

    # Get database connection
    connection = get_database_connection()

    try:
        # Query database
        new_df = query_database(connection)

        # Merge data
        final_df = merge_data(existing_df, new_df)

        # Write to file
        final_df.to_csv(TSV_FILE, sep="\t", index=False)

        print(f"Successfully updated {TSV_FILE}")
        print(f"Total records: {len(final_df)}")
        print(f"New records added: {len(final_df) - len(existing_df)}")

    finally:
        connection.close()


if __name__ == "__main__":
    main()
