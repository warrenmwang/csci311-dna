import os
from typing import Callable
import numpy as np
import string
import sys
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
    m = len(s)
    n = len(t) 
    ret_arr = np.zeros([n+1, m+1])
    result = 0
    for j in range(1,m+1):
        for i in range(1,n+1):
            if s[j-1] == t[i-1]:
                ret_arr[i][j] = ret_arr[i-1][j-1]+ 1
                result = max(result, ret_arr[i][j])
            else:
                ret_arr[i][j] = 0
    return int(result)


def LongestCommonSubsequence(s : str, t : str) -> str:
    
    m = len(s) #j
    n = len(t) #i
    ret_arr = np.zeros([n+1, m+1])
  
    for j in range(1,m+1):
        for i in range(1,n+1):
            if s[j-1] == t[i-1]:
                ret_arr[i][j] = ret_arr[i-1][j-1] + 1
            else:
                ret_arr[i][j] = max(ret_arr[i-1][j],ret_arr[i][j-1])

    return int(ret_arr[n][m])

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
def checkSequences(s : str):
    hold = list(string.ascii_uppercase)
    hold.remove("A")
    hold.remove("T")
    hold.remove("C")
    hold.remove("G")
    if any(x in s for x in hold):
        return False
    else:
        return True
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
                sequence = line.strip().upper()
                if(checkSequences(sequence) == True):
                    D[name] += sequence
                else:
                    print("This is not a valid sequence\nInvalid Sequence Name: " + name)
    return D

def readQuerySequence(filename : str) -> str:
    with open(filename, "r") as f:
        query = f.read()
        query = query.replace("\n", "") # remove newlines
        query.upper()
        if(checkSequences(query)):
            return query
        else:
            print("The provided query contained invalid characters")
            return 

if __name__ == "__main__":
    similarityAlgorithms = {
        "LongestCommonSubstring": LongestCommonSubstring,
        "LongestCommonSubsequence": LongestCommonSubsequence,
        "EditDistance" : EditDistance,
        "NeedleMan_Wunsch" : NeedleMan_Wunsch
    }
    #This gets the dictionary of the sequences
    all_sequences_filename = input("Name of file for database of sequences (in cwd, and in FASTA format): ")
    while not os.path.exists(all_sequences_filename):
        all_sequences_filename = input("File not found. Try again: ")
    D = readDNASequences(all_sequences_filename)
    #This gets the compare sequence
    query_filename = input("Name of query file: ")
    while not os.path.exists(query_filename):
        query_filename = input("File not found. Try again: ")
    query_sequence = readQuerySequence(query_filename)
    
    flag = True
    while(flag == True):
        print("What algorithm would you like to use to compute the similarities between your unknown sequence and those provided?")
        print("".join([f"{i}. {s}\n" for i,s in enumerate(similarityAlgorithms.keys())])) #This just prints the dictionary, why is this so complicated 



        algo_choice = int(input("Select the number: "))
        while(algo_choice not in [0,1,2,3]):
            algo_choice = int(input("Select the number: "))
        if(algo_choice == 0):
            curr_longest = [-1, None]
            for i in D:
                hold = LongestCommonSubstring(D[i], query_sequence)
                print("Comparing String 1: " + i + " to String 2: Query Sequence")
                print("Longest Common Substring: ", hold)
                if(hold > curr_longest[0]):
                    curr_longest[0] = hold
                    curr_longest[1] = i
            print("The longest common substring for the query sequence was: " + curr_longest[1] + " at a total of ", curr_longest[0], "characters")
            #LongestCommonSubstring
            pass
        elif(algo_choice == 1):
            curr_longest = [-1, None]
            for i in D:
                hold = LongestCommonSubsequence(D[i], query_sequence)
                print("Comparing String 1: " + i + " to String 2: Query Sequence")
                print("Longest Common Subsequence: ", hold)
                if(hold > curr_longest[0]):
                    curr_longest[0] = hold
                    curr_longest[1] = i
            print("The longest common sequence for the query sequence was: " + curr_longest[1] + " at a total of ", curr_longest[0], "characters")
        elif(algo_choice == 2):
            curr_longest = [pow(2,31), None]
            for i in D:
                hold = EditDistance(D[i], query_sequence)
                print("Comparing Edit Distance for String 1: " + i + " to String 2: Query Sequence")
                print("Edit Distance is: ", hold)
                print(curr_longest[0])
                if(hold < curr_longest[0]):
                    
                    curr_longest[0] = hold
                    curr_longest[1] = i
            print("The lowest edit distance for the query sequence was: " + curr_longest[1] + " at a total of", curr_longest[0], "characters")
            #LongestCommonSubstring
            #EditDistance
            pass
        elif(algo_choice == 3):
            #NeedleMan_Wunsch
            pass
        val = input("Would you like to do another operation? Type Y to repeat or N to exit")
        print(val)
        if val.strip() == 'n' or val.strip() == 'N':
            break
   
    
    #compute the most similar sequence to query
    #mostSimSeq, sim = FindMostSimilarSeq(query_sequence, D, algo)
    #print(f"The most similar sequence is: {mostSimSeq}\nSimilarity: {sim}")    
