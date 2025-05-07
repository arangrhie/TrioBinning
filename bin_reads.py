import gzip
import sys
import os

def bin_sequences(fastq_file, assignment_file, cutoff):
    bin_files = { }
    bin = ""
    lineNum = 0;
    
    with gzip.open(fastq_file, 'rt') as fq, open(assignment_file, 'r') as af:
        while True:
            fqLine = fq.readline().strip()
            if not fqLine:
                print("Reached the end of fastq file.", file=sys.stderr)
                break

            # Check if the fq line is a header
            # There are qual lines beginning with @, so we need to check the line number
            if fqLine[0] == '@' and lineNum % 4 == 0:
                fqName  = fqLine[1:].split()[0]
                lineNum = 0
            
                # Read the assignment file
                afLine = af.readline().strip()
                if not afLine:
                    print("Reached the end of assignment file.", file=sys.stderr)
                    sys.exit(1)
                    break
                
                parts = afLine.split()
                seqName = parts[0]
                bin = parts[-1]

                if fqName != seqName:
					
                    print(f"Error: sequence names do not match. {bin} :: fq name {fqName} != assignment_file {seqName}", file=sys.stderr)
                    sys.exit(1)
            
                if bin not in bin_files:
                    bin_files[bin] = gzip.open(f'{os.path.splitext(fastq_file)[0]}.{bin}.{cutoff}.fq.gz', 'wt')
            
            bin_files[bin].write(f"{fqLine}\n")
            lineNum += 1
    
    for f in bin_files.values():
        f.close()

def main():
    if len(sys.argv) < 4:
        print("Usage: python bin_reads.py <fastq> <assigned> <cutoff>")
        sys.exit(1)
    
    fastq_file = sys.argv[1]
    assignment_file = sys.argv[2]
    cutoff = float(sys.argv[3])
    
    bin_sequences(fastq_file, assignment_file, cutoff)


if __name__ == "__main__":
    main()