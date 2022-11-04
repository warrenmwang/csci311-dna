import os
import numpy as np
import string
import sys 


def LongestCommonSubstring(s : str, t : str):
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

def LongestCommonSubsequence(s : str, t : str):
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
    m = len(s) + 1
    n = len(t) + 1
    
    d = [[0 for j in range(n)] for i in range(m)]
    
    for i in range(1, m):
        d[i][0] = i
        
    for j in range(1, n):
        d[0][j] = j
        
    for j in range(1, n):
        for i in range(1, m):
            if s[i-1] == t[j-1]:
                substitutionCost = 0
            else:
                substitutionCost = 1

            d[i][j] = min(d[i-1][j]+ 1,            # deletion
                              d[i][j-1] + 1,                # insertion
                              d[i-1][j-1] + substitutionCost)              # substitution
            
    return d[m-1][n-1]


def Needleman_Wunsch(s : str, t : str):
    # based on this pseudocode -> https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm
    EQUAL = 1
    DIFFERENT = -1
    GAP_PENALTY = -5

    m = len(s)
    n = len(t)
   
    d = np.zeros([m+1, n+1])
    for i in range(n):
        d[0][i] = d[0][i-1] + GAP_PENALTY
    for j in range(m):
        d[j][0] = d[j-1][0] + GAP_PENALTY

    for i in range(m + 1):
        for j in range(n + 1):
            maxpoint = -1 * sys.maxsize
            
            if i == 0 and j == 0:
                d[i][j] = 0
            else: 
                if i != 0 and j != 0:
                    if s[i-1] == t[j-1]:
                        matchpoint = EQUAL
                    else:
                        matchpoint = DIFFERENT
                    maxpoint = max(d[i-1][j-1] + matchpoint, maxpoint)
                if i != 0:
                    maxpoint = max(d[i-1][j] + GAP_PENALTY, maxpoint)
                if j != 0:
                    maxpoint = max(d[i][j-1] + GAP_PENALTY, maxpoint)
                d[i][j] = maxpoint
    return int(d[m][n])


def SmithWaterman(s : str, t : str, match = 2, mismatch = -2, gap = -1) -> int:
    dp_matrix = np.zeros((len(s)+1, len(t)+1))
    
    max_val = -1

    for i in range(1, len(s)):
        for j in range(1, len(t)):

            vertical_value = dp_matrix[i-1,j] + gap

            horizontal_value = dp_matrix[i,j-1] + gap

            # s[i-1] is checked with t[j-1] because our indexing in the dp_matrix starts at 1 but the actual string indexing starts at 0. 
            if  s[i-1] == t[j-1]:  
                diagonal_value = dp_matrix[i-1,j-1] + match
            else:
                diagonal_value = dp_matrix[i-1,j-1] + mismatch    

            # if score is negative, replace it with a zero.     
            dp_matrix[i,j] = max(0, diagonal_value, vertical_value, horizontal_value)

            max_val = max(max_val, dp_matrix[i,j])
    return int(max_val)


def checkSequences(s : str):
    hold = list(string.printable)
    hold.remove("A")
    hold.remove("T")
    hold.remove("C")
    hold.remove("G")
    hold.remove("a")
    hold.remove("t")
    hold.remove("c")
    hold.remove("g")
    if any(x in s for x in hold):
        return False
    else:
        return True

def readDNASequences(filename : str) -> dict:
    with open(filename, "r") as f:
        # parse through text and create dict D, (k,v) = (name, sequence)
        lines = f.readlines()
        D = {}
        if lines == "" or lines == []:
            return D
        if lines[0][0] != ">":
            return D
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
        query = query.upper()
        if(query == ""):
            return None
        query = query.replace("\n", "") # remove newlines
        if(checkSequences(query)):
            return query
        else:
            print("The provided query contained invalid characters")
            return None

if __name__ == "__main__":
    LINEBREAK = "-------------------------------------------------------------\n"
    similarityAlgorithms = ["LongestCommonSubstring","LongestCommonSubsequence","EditDistance","NeedleMan_Wunsch","SmithWaterman"]
    
    # user inputs sequences and query file
    all_sequences_filename = input("Name of file for database of sequences (in cwd, and in FASTA format): ")
    D = None
    while D == None or D == {}:
        while not os.path.exists(all_sequences_filename):
            all_sequences_filename = input("File not found. Try again: ")
        D = readDNASequences(all_sequences_filename)
        if(D != {}):
            break
        all_sequences_filename = input("File not found. Try again: ")    
    
    query_filename = input("Name of query file: ")
    query_sequence = None
    while query_sequence == None:
        while not os.path.exists(query_filename):
            query_filename = input("File not found. Try again: ")
        query_sequence = readQuerySequence(query_filename)
        if(query_sequence != None):
            break
        query_filename= input("Provide a file with a valid query sequence: ")
    
    # user selects algo to run on input
    while(True):
        print("What algorithm would you like to use to compute the similarities between your unknown sequence and those provided?")
        print("".join([f"{i}. {s}\n" for i,s in enumerate(similarityAlgorithms)]))
        algo_choice = input("Select the number corresponding to the algorithm you want to run: ")
        listOfAlgoChoices = [str(i) for i in range(len(similarityAlgorithms))]
        while(algo_choice not in listOfAlgoChoices):
            algo_choice =   input("Select the number corresponding to the algorithm you want to run: ")
        algo_choice = int(algo_choice)
        # LongestCommonSubstring
        if(algo_choice == 0):
            curr_longest = [-1, None]
            for i in D:
                hold = LongestCommonSubstring(D[i], query_sequence)
                print("Comparing String 1: " + i + " to String 2: Query Sequence")
                print("Longest Common Substring: ", hold)
                print(LINEBREAK)
                if(hold > curr_longest[0]):
                    curr_longest[0] = hold
                    curr_longest[1] = i
            print("The best match using longest common substring for the query sequence was: " + curr_longest[1] + " at a total of ", curr_longest[0], "characters")
            print(LINEBREAK)
        # LongestCommonSubsequence
        elif(algo_choice == 1):
            curr_longest = [-1, None]
            for i in D:
                hold = LongestCommonSubsequence(D[i], query_sequence)
                print("Comparing String 1: " + i + " to String 2: Query Sequence")
                print("Longest Common Subsequence: ", hold)
                print(LINEBREAK)
                if(hold > curr_longest[0]):
                    curr_longest[0] = hold
                    curr_longest[1] = i
            print("The best match using longest common sequence for the query sequence was: " + curr_longest[1] + " at a total of ", curr_longest[0], "characters")
            print(LINEBREAK)
        # Edit Distance
        elif(algo_choice == 2):
            curr_longest = [pow(2,31), None]
            for i in D:
                hold = EditDistance(D[i], query_sequence)
                print("Comparing String 1: " + i + " to String 2: Query Sequence")
                print("Edit Distance is: ", hold)
                print(LINEBREAK)
                if(hold < curr_longest[0]):
                    curr_longest[0] = hold
                    curr_longest[1] = i
            print("The best match using edit distance for the query sequence was: " + curr_longest[1] + " with the lowest number of ", curr_longest[0], "edits")
            print(LINEBREAK)
        # Needleman-Wunsch
        elif(algo_choice == 3):
            curr_longest = [pow(-2,31), None]
            for i in D:
                hold = Needleman_Wunsch(D[i], query_sequence)
                print("Comparing String 1: " + i + " to String 2: Query Sequence")
                print("Difference : ", hold)
                print(LINEBREAK)
                if(hold > curr_longest[0]):
                    curr_longest[0] = hold
                    curr_longest[1] = i
            print("The best match using Needleman-Wunsch for the query sequence was: " + curr_longest[1] + " with the best score of ", curr_longest[0])
            print(LINEBREAK)
        # Smith-Waterman
        elif(algo_choice == 4):
            curr_longest = [pow(-2,31), None]
            for i in D:
                hold = SmithWaterman(D[i], query_sequence)
                print("Comparing String 1: " + i + " to String 2: Query Sequence")
                print("Score : ", hold)
                print(LINEBREAK)
                if(hold > curr_longest[0]):
                    curr_longest[0] = hold
                    curr_longest[1] = i
            print("The best match using Smith-Waterman for the query sequence was: " + curr_longest[1] + " with the best score of ", curr_longest[0])
            print(LINEBREAK)

        val = input("Would you like to do another operation? Type Y to repeat or N to exit: ")
        print(LINEBREAK)
        if val.strip() == 'n' or val.strip() == 'N':
            break
