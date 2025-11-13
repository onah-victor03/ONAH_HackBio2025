"""
HackBio Internship - Stage 1 Python Task
Team: Glycine

# Task: Write a Python function for calculating the hamming distance between your slack username ('Onah Victor') and twitter/X handle (@onahvictor606@)
# Author: Onah Victor
"""

def hamming_distance(str1, str2):
    # pad shorter string with spaces
    max_len = max(len(str1), len(str2))
    str1 = str1.ljust(max_len)
    str2 = str2.ljust(max_len)

    distance = sum(1 for a, b in zip(str1, str2) if a != b)
    return distance

if __name__ == "__main__":
    slack = "Onah Victor"
    twitter = "@onahvictor606@"
    print("Slack:", slack)
    print("Twitter:", twitter)
    print("Hamming distance:", hamming_distance(slack, twitter))
