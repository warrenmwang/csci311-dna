import os
from typing import Callable

def FindMostSimilarSeq(t : str, D : dict, algo : Callable) -> tuple[str, int]:
    bestSim = -10e-10
    bestSeq = ""
    for i in range(len(D)):
        sim = ComputeSimilarity(t, D[i], algo)
        if sim > bestSim:
            bestSim = sim
            bestSeq = D[i]
    return bestSeq, bestSim

def ComputeSimilarity(t : str, s : str, algo : Callable) -> int:
    # Compute the similarity between t and s
    # using the function passed in as algo
    # TODO
    pass

def LongestCommonSubstring(s : str, t : str) -> str:
    # Compute the longest common substring between s and t
    # TODO
    pass

def LongestCommonSubsequence(s : str, t : str) -> str:
    # Compute the longest common subsequence between s and t
    # TODO
    pass

def EditDistance(s : str, t : str) -> int:
    # Compute the Levenshtein distance between s and t
    '''
    created from pseudocode found here: https://en.wikipedia.org/wiki/Levenshtein_distance
    '''
    m = len(s)
    n = len(t)
    
    d = [[0 for j in range(n)] for i in range(m)]
    
    for i in range(m):
        d[i][0] = i
        
    for j in range(n):
        d[0][j] = j
        
    for j in range(n):
        for i in range(m):
            if s[i] == t[j]:
                substitutionCost = 0
            else:
                substitutionCost = 1
                
            d[i][j] = min(d[i-1][j]   + 1,                # deletion
                          d[i][j-1]   + 1,                # insertion
                          d[i-1][j-1] + substitutionCost) # substitution
            
    return d[m-1][n-1]

def NeedleMan_Wunsch(s : str, t : str) -> int:
    # TODO
    pass

def readDNASequences(filename : str) -> dict:
    with open(filename, "r") as f:
        # parse through text and create dict D, (k,v) = (name, sequence)
        lines = f.readlines()
        D = {}
        for line in lines:
            if line[0] == ">":
                name = line[1:].strip()
                D[name] = ""
            else:
                D[name] += line.strip()
    return D

def readQuerySequence(filename : str) -> str:
    with open(filename, "r") as f:
        query = f.read()
        query = query.replace("\n", "") # remove newlines
    return query

if __name__ == "__main__":
    similarityAlgorithms = {
        "LongestCommonSubstring": LongestCommonSubstring,
        "LongestCommonSubsequence": LongestCommonSubsequence,
        "EditDistance" : EditDistance,
        "NeedleMan_Wunsch" : NeedleMan_Wunsch
    }

    all_sequences_filename = input("Name of file for database of sequences (in cwd, and in FASTA format): ")
    while not os.path.exists(all_sequences_filename):
        all_sequences_filename = input("File not found. Try again: ")
    D = readDNASequences(all_sequences_filename)

    query_filename = input("Name of query file: ")
    while not os.path.exists(query_filename):
        query_filename = input("File not found. Try again: ")
    query_sequence = readQuerySequence(query_filename)

    print("What algorithm would you like to use to compute the similarities between your unknown sequence and those provided?")
    print("".join([f"{i}. {s}\n" for i,s in enumerate(similarityAlgorithms.keys())]))
    algo_choice = int(input("Select the number: "))
    algo = similarityAlgorithms[list(similarityAlgorithms.keys())[algo_choice]]

    #compute the most similar sequence to query
    mostSimSeq, sim = FindMostSimilarSeq(query_sequence, D, algo)
    print(f"The most similar sequence is: {mostSimSeq}\nSimilarity: {sim}")