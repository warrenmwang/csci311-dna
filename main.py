def FindMostSimilarSeq(t : str, D : dict) -> tuple[str, int]:
    bestSim = -10e-10
    bestSeq = ""
    for i in range(len(D)):
        sim = ComputeSimilarity(t, D[i])
        if sim > bestSim:
            bestSim = sim
            bestSeq = D[i]
    return bestSeq, bestSim

def ComputeSimilarity(t : str, s : str) -> int:
    # Compute the similarity between t and s
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
    # Compute the edit distance between s and t
    # number of edit operations required to transform
    # s into t
    # TODO
    pass

def NeedleMan_Wunsch(s : str, t : str) -> int:
    # TODO
    pass

if __name__ == "__main__":
    with open("DNA_sequences.txt", "r") as f:
        D = f.readlines()
        # TODO: parse through text and create dict D
        # format to store, TBD

    with open("DNA_query.txt", "r") as f:
        query = f.read()
        query = query.replace("\n", "") # remove newlines

    # compute the most similar sequence to query
    mostSimSeq, sim = FindMostSimilarSeq(query, D)
    print(mostSimSeq, sim)