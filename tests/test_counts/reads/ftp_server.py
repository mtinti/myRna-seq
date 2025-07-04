#!/usr/bin/env python3

from pyftpdlib.authorizers import DummyAuthorizer
from pyftpdlib.handlers import FTPHandler
from pyftpdlib.servers import FTPServer
import os
import glob

# Define port number
PORT = 2121

def main():
    # Use current directory
    directory = os.getcwd()
    
    # Create a dummy authorizer for managing users
    authorizer = DummyAuthorizer()
    
    # Add anonymous user with read permissions
    authorizer.add_anonymous(directory, perm='elr')  # Read-only permissions

    # Create handler
    handler = FTPHandler
    handler.authorizer = authorizer
    
    # Create server
    server = FTPServer(("0.0.0.0", PORT), handler)
    
    # Get hostname
    hostname = "localhost"
    
    # Display server info
    print(f"Anonymous FTP server starting on port {PORT}")
    print(f"Serving files from current directory: {directory}")
    
    # Find FASTQ files in the current directory
    fastq_files = glob.glob("*.fastq.gz") + glob.glob("*/*.fastq.gz")
    print(fastq_files)
    if fastq_files:
        print("\nAvailable FASTQ files:")
        for file in fastq_files:
            print(f"- ftp://{hostname}:{PORT}/{file}")
        
        # Sample CSV entries
        print("\nSample CSV entries for these files:")
        for i, file in enumerate(fastq_files):
            sample_name = f"SAMPLE_{i+1:02d}"
            if file.endswith(".1.fastq.gz") and file.replace(".1.fastq.gz", ".2.fastq.gz") in fastq_files:
                # Paired-end file
                r1 = file
                r2 = file.replace(".1.fastq.gz", ".2.fastq.gz")
                print(f"{sample_name},paired,ftp,ftp://{hostname}:{PORT}/{r1},ftp://{hostname}:{PORT}/{r2},,")
            elif file.endswith(".2.fastq.gz"):
                # This is a R2 file, already handled with its R1 pair
                continue
            else:
                # Single-end file
                print(f"{sample_name},single,ftp,ftp://{hostname}:{PORT}/{file},,,")
    else:
        print("\nNo FASTQ files found in the current directory.")
    
    print("\nPress Ctrl+C to stop the server")
    
    # Start the server
    server.serve_forever()

if __name__ == "__main__":
    main()