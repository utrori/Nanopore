failed = 0
passed = 0
failed_files = 0
passed_files = 0
with open('fast5.listing.txt') as f:
    for line in f:
        row = line.split()
        size = int(row[2])
        if 'fast5_pass' in line:
            passed += size
            passed_files += 1
        elif 'fast5_fail' in line:
            failed_files += 1
            failed += size
print(passed)
print(failed)
print(passed_files)
print(failed_files)

