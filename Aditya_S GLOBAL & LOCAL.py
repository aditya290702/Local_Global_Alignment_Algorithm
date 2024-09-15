#GLOBAL ALIGNMENT

import numpy as np
def Global():
    # My scoring schema
    match = 1
    mismatch = -1
    gap = -2

    Seq = "ATGC"
    Qr = "AATGG"
    Seq = list(Seq)
    Qr = list(Qr)
    Seq.insert(0, 0)
    Qr.insert(0, 0)
    S_len = len(Seq)
    Q_len = len(Qr)

    # Matrix Creation
    Matrix = np.zeros([S_len, Q_len], dtype=int)
    for a in range(S_len):
        # print(Seq[a], end="  ")
        for b in range(Q_len):
            Matrix[a, 0] = a * gap
            Matrix[0, b] = b * gap

    # Score Implementation in Matrix
    for a in range(1, S_len):
        for b in range(1, Q_len):
            if Seq[a] == Qr[b]:
                row_score = Matrix[a - 1, b] + gap
                col_score = Matrix[a, b - 1] + gap
                diagonal_score = Matrix[a - 1, b - 1] + match
                Matrix[a, b] = max(row_score, col_score, diagonal_score)
            else:
                row_score = Matrix[a - 1, b] + gap
                col_score = Matrix[a, b - 1] + gap
                diagonal_score = Matrix[a - 1, b - 1] + mismatch
                Matrix[a, b] = max(row_score, col_score, diagonal_score)

    # To Print the Matrix
    print("     ", end="")
    for ele in Qr:
        print(ele, end="    ")
    print()
    for a in range(S_len):
        print(Seq[a], end="  ")
        for b in range(Q_len):
            Matrix[a, 0] = a * gap
            Matrix[0, b] = b * gap
            print(f"{Matrix[a, b]: >3}", end="  ")
        print()

    # Tracing
    Aligned = []
    Query = []
    Seq_iter = S_len - 1
    Query_iter = Q_len - 1

    while Seq_iter > 0 and Query_iter > 0:
        current_score = Matrix[Seq_iter, Query_iter]
        diagonal_score = Matrix[Seq_iter - 1, Query_iter - 1]
        left_score = Matrix[Seq_iter, Query_iter - 1]
        up_score = Matrix[Seq_iter - 1, Query_iter]

        if Seq[Seq_iter] == Qr[Query_iter]:
            score = match
        else:
            score = mismatch

        if current_score == diagonal_score + score:
            Aligned.append(Seq[Seq_iter])
            Query.append(Qr[Query_iter])
            Seq_iter = Seq_iter - 1
            Query_iter = Query_iter - 1
        elif current_score == up_score + gap:
            Aligned.append(Seq[Seq_iter])
            Query.append("_")
            Seq_iter = Seq_iter - 1
        elif current_score == left_score + gap:
            Aligned.append("_")
            Query.append(Qr[Query_iter])
            Query_iter = Query_iter - 1

    # Handle remaining elements if any
    while Seq_iter > 0:
        Aligned.append(Seq[Seq_iter])
        Query.append("_")
        Seq_iter -= 1

    while Query_iter > 0:
        Aligned.append("_")
        Query.append(Qr[Query_iter])
        Query_iter -= 1

    # Print the aligned sequences
    print("Aligned Sequence:", ''.join(Aligned[::-1]))
    print("Aligned Query:   ", ''.join(Query[::-1]))


#-------------------------------------------------------------------------------------

def local():
    import numpy as np
    # LOCAL ALIGNMENT

    S = "ATGCGCTTGC"
    Q = "AGCTG"
    Sequence = "_" + S
    Query = "_" + Q
    lens1 = len(Sequence)
    lens2 = len(Query)
    L = np.zeros((lens1, lens2), dtype=int)
    Alignment = 0

    print("  ", end="")
    for ele in Query:
        print(ele, end=" ")
    print()

    for i in range(1, lens1):
        for j in range(1, lens2):
            if Sequence[i] == Query[j]:
                match = 1
                L[i, j] = match
            else:
                L[i, j] = 0

    for i in range(lens1):
        print(Sequence[i], end=" ")
        for j in range(lens2):
            print(L[i, j], end=" ")
        print()

    seq_list = []
    query_list = []

    while len(seq_list) < len(S):
        if L[i, j] == 1:
            seq_list.append(Sequence[i])
            query_list.append(Query[j])
            i = i - 1
            j = j - 1
        elif L[i, j] == 0 and L[i - 1, j] == 1 or L[i, j - 1] == 1:
            seq_list.append(Sequence[i])
            query_list.append("_")
            i = i - 1
        elif L[i, j] == 0 and L[i - 1, j] == 0 or L[i, j - 1] == 0:
            seq_list.append(Sequence[i])
            query_list.append("_")
            i = i - 1

    # Print the reversed lists
    print(seq_list[::-1])
    print(query_list[::-1])


def main():
    print("GlobaL ALIGNMENT")
    Global()
    print()
    print()
    print()
    print("LocaL ALIGNMENT")
    local()


if __name__ == '__main__':
    main()
