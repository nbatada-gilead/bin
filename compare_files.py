#!/usr/bin/env python3
import os
import filecmp
import difflib
import sys
import argparse

def compare_files(master_file, slave_file):
    if not os.path.exists(master_file):
        return f"Error: {master_file} does not exist", None
    if not os.path.exists(slave_file):
        return f"Error: {slave_file} does not exist", None

    with open(master_file, 'r') as f1, open(slave_file, 'r') as f2:
        diff = list(difflib.unified_diff(f1.readlines(), f2.readlines(), lineterm=''))
        if diff:
            return "Different", diff
        else:
            return "Identical", None

def generate_report(master_dir, slave_dir):
    if not os.path.exists(master_dir):
        print(f"Error: {master_dir} does not exist")
        return
    if not os.path.exists(slave_dir):
        print(f"Error: {slave_dir} does not exist")
        return

    report = []
    master_files = {f for f in os.listdir(master_dir) if f.endswith(('.py', '.txt', '.sh'))}
    slave_files = {f for f in os.listdir(slave_dir) if f.endswith(('.py', '.txt', '.sh'))}

    all_files = master_files.union(slave_files)

    for file_name in all_files:
        master_file = os.path.join(master_dir, file_name)
        slave_file = os.path.join(slave_dir, file_name)

        if file_name not in master_files:
            report.append(f"{file_name}\tMissing in Master")
        elif file_name not in slave_files:
            report.append(f"{file_name}\tMissing in Slave")
        else:
            status, diff = compare_files(master_file, slave_file)
            if status == "Identical":
                report.append(f"{file_name}\tIdentical")
            else:
                diff_summary = []
                for line in diff:
                    if line.startswith('+') and not line.startswith('+++'):
                        diff_summary.append(f"Master: {line.strip()}")
                    elif line.startswith('-') and not line.startswith('---'):
                        diff_summary.append(f"Slave: {line.strip()}")
                
                if diff_summary:
                    report.append(f"{file_name}\tDifferent\t{', '.join(diff_summary[:3])}")

    for line in report:
        print(line)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare files between master and slave directories.")
    parser.add_argument('-m', '--master', required=True, help="Path to the master directory")
    parser.add_argument('-s', '--slave', required=True, help="Path to the slave directory")
    args = parser.parse_args()

    generate_report(args.master, args.slave)
