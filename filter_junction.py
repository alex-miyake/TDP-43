import pyarrow.parquet as pq
import pyarrow.compute as pc
import pyarrow.csv as pv
import pyarrow as pa
import pandas as pd
import zipfile
import os

# Extract all files from the ZIP
def extract_all_parquet_from_zip(zip_path, extract_dir='extracted_parquet'):
    """
    Extracts all Parquet files from a ZIP archive.

    Parameters:
    zip_path (str): Path to the ZIP file containing parquet files.
    extract_dir (str): Directory where files will be extracted (default: 'extracted_parquet').

    Returns:
    List[str]: List of paths to the extracted Parquet files.

    Raises:
    Exception: If no parquet files are found or if extraction fails.
    """
    
    try:
        print(f"Extracting files from {zip_path}...")
        
        # Create directory to extract files
        if not os.path.exists(extract_dir):
            os.makedirs(extract_dir)
        
        # Extract all parquet files
        with zipfile.ZipFile(zip_path, 'r') as z:
            z.extractall(extract_dir)
        
        # Get list of all parquet files
        parquet_files = [os.path.join(extract_dir, f) for f in os.listdir(extract_dir) if f.endswith('.parquet')]
        
        if len(parquet_files) == 0:
            raise Exception(f"No parquet files found in the extracted directory {extract_dir}.")
        
        print(f"Successfully extracted {len(parquet_files)} parquet files.")
        return parquet_files
    
    except zipfile.BadZipFile:
        print("Error: Failed to unzip the file. It might be corrupted.")
        raise
    except Exception as e:
        print(f"Error during extraction: {e}")
        raise

# Load Parquet files
def load_parquet_files(parquet_files):
    """
    Loads Parquet files and returns them as pyarrow.Table.

    Parameters:
    parquet_files (List[str]): List of paths to the Parquet files.

    Returns:
    pyarrow.Table: Loaded parquet table.

    Raises:
    Exception: If there is an issue reading the parquet files.
    """
    try:
        print("Loading parquet files...")
        
        for parquet_file in parquet_files:
            print(f"Processing file: {parquet_file}")
            parquet_file = pq.ParquetFile(parquet_file, thrift_string_size_limit = 1000000000)
            parquet_dataset = parquet_file.read()
            print("Successfully read parquet file")
            print(parquet_file.metadata)
            print(parquet_file.schema)
        print("Successfully loaded parquet files.")
        
        return parquet_dataset

    except Exception as e:
        print(f"Error during parquet loading: {e}")
        raise

# Filter splice events function
def filter_splice(splice_events_path):
    """
    Filters the splice events DataFrame based on specific cryptic event categories.

    Parameters:
    splice_events_path (str): Path to the CSV file containing splice events.

    Returns:
    pandas.DataFrame: Filtered DataFrame with cryptic splice events.

    Raises:
    Exception: If there is an issue during filtering.
    """
    # list of cryptic junction categories
    options = ['novel_acceptor','novel_exon_skip','novel_donor']
    # filter based on this list
    try:
        splice_df = pd.read_csv(splice_events_path)
        filtered_df = splice_df.loc[splice_df['junc_cat'].isin(options)]
        print("successfully filtered splice events to cryptic events")
        return filtered_df
    
    except Exception as e:
        print(f"Error during splice filtering: {e}")
        raise

# Join datasets together 
def join_datasets(parquet_table, splice_events, metadata):
    """
    Converts splice events pandas DataFrame to PyArrow table.
    Joins the Parquet table with the splice events and metadata.

    Parameters:
    parquet_table (pyarrow.Table): Pyarrow Table of the Parquet dataset.
    splice_events (pandas.DataFrame): DataFrame of filtered splice events.
    metadata: CSV file containing sample metadata.

    Returns:
    pyarrow.Table: Joined table containing data from both datasets.

    Raises:
    Exception: If there is an issue during the join process.
    """
    try:
        ds2 = pa.Table.from_pandas(splice_events)
        splice_joined_ds = parquet_table.join(ds2, keys="junction_coords", join_type="inner")
        print("Successfully joined splice events to hek junctions")

    except Exception as e:
        print(f"Error joining splice events with hek junctions: {e}")
        raise

    try:
        ds3 = pv.read_csv(metadata)
        all_added_ds = splice_joined_ds.join(ds3, keys="sample_name")      
        print("Successfully joined metadata to hek junctions")

    except Exception as e:
        print(f"Error joining  with hek junctions: {e}")
        raise

    return all_added_ds



def count_joined(joined_df):
    """    
    Creates a summarisation table with the number of cryptics per genotype, along with number of crpytics associated with TDP43 knockdown, and number of events where rescue is induced.
    
    Parameters:
    joined_df (DataFrame)

    Returns:
    counted_df (DataFrame): table with genotype, cryptic count, TDP43 count, and resue count.

    Raises:
    Exception: If there is an issue during counting.
    """
    try:    
        joined_df = joined_df.to_pandas()
        counted_df = joined_df.groupby('genotype').agg(cryptic_count=('genotype', 'size'),TDP34_kd_associated=('TDP43_kd', lambda x: (x == 'siTDP43').sum()),rescueInduced=('rescueExpression', lambda x: (x == 'rescueInduced').sum()))
        print('Successfully counted number of cryptics per gene')
        return counted_df
    
    except Exception as e:
        print(f"Error counting joined junctions: {e}")
    raise


# filter 
def filter_joined(joined_df):
    """
    Filters the joined dataset for only entries with TDP43 knockdown, and if rescue is induced.

    Parameters:
    joined_df (pyarrow.Table): Joined table containing both HEK junctions and splice events.

    Returns:
    pyarrow.Table: Filtered table.

    Raises:
    Exception: If there is an issue during filtering.
    """
    # Filter for only entries with TDP43 knockdown
    try:
        filtered_df = joined_df.filter(pc.equal(joined_df["TDP43_kd"], "siTDP43")) 
        print(f"Result: There are {filtered_df.num_rows} rows of cryptic events that appear with TDP43_kd")
        pv.write_csv(filtered_df, "unfiltered_cryptic_data.csv")

    except Exception as e:
        print(f"Error filtering with joined junctions: {e}")
        raise

    # filter for only entries with RescueInduced
    try:
        filtered2_df = filtered_df.filter(pc.equal(filtered_df["rescueExpression"], "rescueInduced")) 
        print(f"Result: There are {filtered2_df.num_rows} rows of rescued cryptic events that appear with TDP43_kd")

    except Exception as e:
        print(f"Error filtering with joined junctions (2): {e}")
        raise

    return filtered2_df


# Main processing function
def process_data(zip_path, splice_events, sample_metadata, extract_dir='extracted_parquet'):
    try:
        print("Starting data processing...")
        
        # Extract all Parquet files from the ZIP archive
        parquet_files = extract_all_parquet_from_zip(zip_path, extract_dir=extract_dir)
        parquet_dataset = load_parquet_files(parquet_files)

        # Filter splice events csv
        filtered_df = filter_splice(splice_events)
        joined_df = join_datasets(parquet_dataset, filtered_df, sample_metadata)
        filtered_data_df = filter_joined(joined_df)
        counted_data_df = count_joined(joined_df)

        # Clean up by removing the extracted directory
        for file in parquet_files:
            os.remove(file)
        os.rmdir(extract_dir)
        print("Temporary files cleaned up successfully.")
        
        return filtered_data_df, counted_data_df
    
    except Exception as e:
        print(f"Error during data processing: {e}")
        raise


zip_path = 'input_data/hek_all_junctions.parquet.zip'  # Parquet ZIP file
splice_events_path = 'input_data/splice_events.csv'  # Splice events CSV file
sample_metadata_path = 'input_data/metadata_halleger_hek.csv'  # Sample metadata CSV file

try:
    # Load Splice Events and Sample Metadata
    print("Loading splice events and sample metadata...")
    
    # Run the processing
    filtered_df, counted_df = process_data(zip_path, splice_events_path, sample_metadata_path)
    pv.write_csv(filtered_df, "filtered_data.csv")
    counted_df.to_csv('cryptic_summarisation.csv')
    print("Results saved successfully.")
    
except Exception as e:
    print(f"Error during the script execution: {e}")