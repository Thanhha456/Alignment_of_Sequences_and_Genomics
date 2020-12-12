"""
We will implement four functions. The first pair of functions will return matrices that we will use in computing the alignment of two sequences.
The second pair of functions will return global and local alignments of two input sequences based on a provided alignment matrix. You will then use these functions in Application 4 to analyze two problems involving comparison of similar sequences.
"""
import urllib.request as urllib2
import random
import matplotlib.pyplot as plt
import math
import numpy

HUMAN_URL = "http://storage.googleapis.com/codeskulptor-alg/alg_HumanEyelessProtein.txt"
FRUITIFLY_URL = "http://storage.googleapis.com/codeskulptor-alg/alg_FruitflyEyelessProtein.txt"
PAM50_URL = "http://storage.googleapis.com/codeskulptor-alg/alg_PAM50.txt"
CONSENSUS_PAX50_URL = "http://storage.googleapis.com/codeskulptor-alg/alg_ConsensusPAXDomain.txt"
CHECKT_LIST = "http://storage.googleapis.com/codeskulptor-assets/assets_scrabble_words3.txt"
scoring_matrix = load_data_table(PAM50_URL)
seq_x = read_file(HUMAN_URL)
seq_y = read_file(FRUITIFLY_URL)
consensus_pax = read_file(CONSENSUS_PAX50_URL)

def build_scoring_matrix(alphabet, diag_score, off_diag_score, dash_score):
    """
    :param alphabet:
    :param diag_score:
    :param off_diag_score:
    :param dash_score:
    :return: dictionary of dictionaries whose entries are indexed by pairs of characters in alphabet plus ’-’
    """
    alphabet_copy = list(alphabet)
    alphabet_copy.append("-")
    scoring_matrix = {}
    for chr_key in alphabet_copy:
        row = {}
        for chr_value in alphabet_copy:
            if chr_value == chr_key and chr_key != "-":
                row[chr_value] = diag_score
            elif chr_value == "-" or chr_key == "-":
                row[chr_value] = dash_score
            else:
                row[chr_value] = off_diag_score
        scoring_matrix[chr_key] = row

    return scoring_matrix

def compute_alphabet(seq_x, seq_y):
    """
    create the alphabet for the two sequences
    """
    alphabet = set()
    for dummy in seq_x:
        alphabet.add(dummy)
    for dummy in seq_y:
        alphabet.add(dummy)
    return list(alphabet)

def compute_alignment_matrix(seq_x, seq_y, scoring_matrix, global_flag):
    """
    Takes as input two sequences seq_x and seq_y whose elements share a common alphabet with the scoring matrix scoring_matrix.
    The function computes and returns the alignment matrix for seq_x and seq_y as described in the Homework. If global_flag is True,
    each entry of the alignment matrix is computed using the method for Algorithm ComputeGlobalAlignmentScores.
    If global_flag is False, each entry is computed using the method for local alignment matrix .
    :param seq_x:
    :param seq_y:
    :param scoring_matrix:
    :param global_flag:
    :return:the alignment matrix for seq_x and seq_y
    """
    alignment_matrix = [[0 for dummy in range(len(seq_y) + 1)] for dummy in range(len(seq_x) + 1)]
    for idx_x in range(1, len(alignment_matrix)):
        alignment_matrix[idx_x][0] = alignment_matrix[idx_x - 1][0] + scoring_matrix[seq_x[idx_x - 1]]["-"]
        if alignment_matrix[idx_x][0] < 0 and not global_flag:
            alignment_matrix[idx_x][0] = 0

    for idx_y in range(1, len(alignment_matrix[0])):
        alignment_matrix[0][idx_y] = alignment_matrix[0][idx_y - 1] + scoring_matrix["-"][seq_y[idx_y - 1]]
        if alignment_matrix[0][idx_y] < 0 and not global_flag:
            alignment_matrix[0][idx_y] = 0

    for idx_xx in range(1, len(alignment_matrix)):
        for idx_yy in range(1, len(alignment_matrix[0])):
            score1 = alignment_matrix[idx_xx - 1][idx_yy - 1] + scoring_matrix[seq_x[idx_xx - 1]][seq_y[idx_yy - 1]]
            score2 = alignment_matrix[idx_xx - 1][idx_yy] + scoring_matrix[seq_x[idx_xx - 1]]["-"]
            score3 = alignment_matrix[idx_xx][idx_yy - 1] + scoring_matrix["-"][seq_y[idx_yy - 1]]
            if max(score1, score2, score3) <= 0 and not global_flag:
                alignment_matrix[idx_xx][idx_yy] = 0
            else:
                alignment_matrix[idx_xx][idx_yy] = max(score1, score2, score3)
    return alignment_matrix

def compute_global_alignment(seq_x, seq_y, scoring_matrix, alignment_matrix):
    """
    Takes as input two sequences seq_x and seq_y whose elements share a common alphabet with the scoring matrix scoring_matrix.
    This function computes a global alignment of seq_x and seq_y using the global alignment matrix alignment_matrix.
    The function returns a tuple of the form (score, align_x, align_y) where score is the score of the global alignment align_x and align_y.
    Note that align_x and align_y should have the same length and may include the padding character '-'.
    :param seq_y:
    :param scoring_matrix:
    :return:score, align_x, align_y
    """
    align_x = ""
    align_y = ""
    idx_i = len((seq_x))
    idx_j = len(seq_y)

    while idx_i != 0 and idx_j != 0:
        if alignment_matrix[idx_i][idx_j] == alignment_matrix[idx_i - 1][idx_j - 1] + scoring_matrix[seq_x[idx_i - 1]][seq_y[idx_j - 1]]:
            align_x = seq_x[idx_i - 1] + align_x
            align_y = seq_y[idx_j - 1] + align_y
            idx_i = idx_i - 1
            idx_j = idx_j - 1
        else:
            if alignment_matrix[idx_i][idx_j] == alignment_matrix[idx_i - 1][idx_j] + scoring_matrix[seq_x[idx_i - 1]]["-"]:
                align_x = seq_x[idx_i - 1] + align_x
                align_y = "-" + align_y
                idx_i = idx_i - 1
            else:
                align_x = "-" + align_x
                align_y = seq_y[idx_j - 1] + align_y
                idx_j = idx_j - 1
    while idx_i != 0:
        align_x = seq_x[idx_i - 1] + align_x
        align_y = "-" + align_y
        idx_i = idx_i - 1
    while idx_j != 0:
        align_x = "-" + align_x
        align_y = seq_y[idx_j - 1] + align_y
        idx_j = idx_j - 1
    return alignment_matrix[-1][-1], align_x, align_y

def compute_local_alignment(seq_x, seq_y, scoring_matrix, alignment_matrix):
    """
    Takes as input two sequences seq_x and seq_y whose elements share a common alphabet with the scoring matrix scoring_matrix.
    This function computes a local alignment of seq_x and seq_y using the local alignment matrix alignment_matrix.
    The function returns a tuple of the form (score, align_x, align_y)
    where score is the score of the optimal local alignment align_x and align_y.
    Note that align_x and align_y should have the same length and may include the padding character '-'
    :param seq_x:
    :param seq_y:
    :param scoring_matrix:
    :param alignment_matrix:
    :return:
    """
    align_x = ""
    align_y = ""
    idx_i = len((seq_x))
    idx_j = len(seq_y)
    max_score = -float("inf")
    for dummy_x in range(len(alignment_matrix)):
        for dummy_y in range(len(alignment_matrix[0])):
            if alignment_matrix[dummy_x][dummy_y] > max_score:
                max_score = alignment_matrix[dummy_x][dummy_y]
                max_x, max_y = dummy_x, dummy_y
    idx_i = max_x
    idx_j = max_y
    while idx_i != 0 and idx_j != 0:
        if alignment_matrix[idx_i][idx_j] == alignment_matrix[idx_i - 1][idx_j - 1] + scoring_matrix[seq_x[idx_i - 1]][seq_y[idx_j - 1]]:
            align_x = seq_x[idx_i - 1] + align_x
            align_y = seq_y[idx_j - 1] + align_y
            # print("align_x1 = ", align_x)
            # print("align_y1 = ", align_y)
            idx_i -= 1
            idx_j -= 1
        else:
            if alignment_matrix[idx_i][idx_j] == alignment_matrix[idx_i - 1][idx_j] + scoring_matrix[seq_x[idx_i - 1]]["-"]:
                align_x = seq_x[idx_i - 1] + align_x
                align_y = "-" + align_y
                # print("align_x2 = ", align_x)
                # print("align_y2 = ", align_y)
                idx_i -= 1
            else:
                align_x = "-" + align_x
                align_y = seq_y[idx_j - 1] + align_y
                # print("align_x3 = ", align_x)
                # print("align_y3 = ", align_y)
                idx_j -= 1

        if alignment_matrix[idx_i][idx_j] == 0:
            break
        elif idx_i == 0 or idx_j == 0:
            break
    count_x = 0
    while len(align_x) != 0:
        if align_x.startswith("-"):
            align_x = align_x[1:]
            count_x += 1
        else:
            break
    count_y = 0
    while len(align_y) != 0:
        if align_y.startswith("-"):
            align_y = align_y[1:]
            count_y += 1
        else:
            break

    if count_x != 0:
        for dummy_x in range(count_x):
            align_y = align_y[1:]
    if count_y != 0:
        for dummy_y in range(count_y):
            align_x = align_x[1:]
    return max_score, align_x, align_y

# #########################################################################################################
def read_word(word_list_url):
    """
    read word data
    """
    data_file = urllib2.urlopen(word_list_url)
    data = data_file.read()
    data = data.decode()
    word_list = data.split('\n')
    return word_list

def read_file(data_url):
    """
    read data file
    """
    data_file = urllib2.urlopen(data_url)
    data = data_file.read()
    data = data.decode()
    data = data.rstrip()
    return data

def load_data_table(data_url):
    """
    Import a table of county-based cancer risk data
    from a csv format file
    """
    data_file = urllib2.urlopen(data_url)
    data = data_file.read()
    data = data.decode()
    data_lines = data.split('\n')
    data_list = []
    for line in data_lines:
        data_list.append(line.split())
    result = {}
    for dummy_x in range(len(data_list[0])):
        row = {}
        for dummy_y in range(len(data_list[0])):
            row[data_list[0][dummy_y]] = int(data_list[dummy_x + 1][dummy_y + 1])
        result[data_list[0][dummy_x]] = row

    return result

def local_align_human_fruitfly(seq_x, seq_y, scoring_matrix):
    """
    compute the local alignments of the sequences of HumanEyelessProtein and FruitflyEyelessProtein using the PAM50 scoring matrix and enter the score and local alignments for these two sequences
    :param seq_x:
    :param seq_y:
    :return: score , hum_ff_local, ff_hum_local
    """
    alignment_matrix = compute_alignment_matrix(seq_x, seq_y, scoring_matrix, global_flag=False)
    #local alignment score of human-eyeless seq and fruitfly seq and their local alignments
    score , hum_ff_local, ff_hum_local = compute_local_alignment(seq_x, seq_y, scoring_matrix, alignment_matrix)
    return score , hum_ff_local, ff_hum_local

def global_align_hum_consus_ff(seq_x, seq_y, consensus_pax, scoring_matrix):
    """
     For each of the two sequences of the local alignment computed in Question 1, do the following:

    Delete any dashes '-' present in the sequence.
    Compute the global alignment of this dash-less sequence with the ConsensusPAXDomain sequence.
    Compare corresponding elements of these two globally-aligned sequences (local vs. consensus)
    and compute the percentage of elements in these two sequences that agree.

    :param seq_x:
    :param seq_y:
    :param consensus_pax:
    :return:
    """
    score, hum_ff_local, ff_hum_local = local_align_human_fruitfly(seq_x, seq_y,scoring_matrix )

    hum_ff_local =hum_ff_local.replace("-", "")
    ff_hum_local = ff_hum_local.replace("-", "")
    alignment_matrix_hum_cons_gl = compute_alignment_matrix(hum_ff_local, consensus_pax, scoring_matrix, global_flag=True)
    alignment_matrix_ff_cons_gl = compute_alignment_matrix(ff_hum_local, consensus_pax, scoring_matrix, global_flag=True)
    hum_cons_gl_score, hum_cons_gl, cons_hum_gl = compute_global_alignment(hum_ff_local, consensus_pax, scoring_matrix, alignment_matrix_hum_cons_gl)
    ff_cons_gl_score, ff_cons_gl, cons_ff_gl = compute_global_alignment(ff_hum_local, consensus_pax, scoring_matrix, alignment_matrix_ff_cons_gl)
    #Compare coresponding elements os these two globally alligned seq( local vs. consensus) and compute the percentage agree elements
    count = 0
    for idx in range(len(hum_cons_gl)):
        if hum_cons_gl[idx] == cons_hum_gl[idx]:
            count += 1
    hum_cons_agree = float(count)/len(hum_cons_gl)
    print((count))
    count = 0
    for idx in range(len(ff_cons_gl)):
        if ff_cons_gl[idx] == cons_ff_gl[idx]:
            count += 1
    print(count)
    ff_cons_agree =float(count)/len(ff_cons_gl)
    # print("hum_cons_agree, ff_cons_agree = ",hum_cons_agree,",", ff_cons_agree)
    return hum_cons_agree, ff_cons_agree

# alignment_matrix = compute_alignment_matrix(seq_x, seq_y, scoring_matrix, global_flag=False)
############################################################################################################

def generate_null_distribution(seq_x, seq_y, scoring_matrix, num_trials):
    """
    This function should return a dictionary scoring_distribution that represents an un-normalized distribution
    generated by performing the following process num_trials times:

    - Generate a random permutation rand_y of the sequence seq_y using random.shuffle().
    - Compute the maximum value score for the local alignment of seq_x and rand_y using the score matrix scoring_matrix.
    - Increment the entry score in the dictionary scoring_distribution by one.""
    Use the function generate_null_distribution to create a distribution with 1000 trials using the protein sequences HumanEyelessProtein and FruitflyEyelessProtein (using the PAM50 scoring matrix).
    Next, create a bar plot of the normalized version of this distribution using plt.bar in matplotlib (or your favorite plotting tool)
    :param seq_x:
    :param seq_y:
    :param scoring_matrix:
    :param num_trials:
    :return:
    """
    scoring_set = {}
    for dummy in range(num_trials):
        ran_y = "".join(random.sample(seq_y, len(seq_y)))
        alignment_matrix = compute_alignment_matrix(seq_x, ran_y, scoring_matrix, global_flag=False)
        score, dummy_x, dummy_y = compute_local_alignment(seq_x, ran_y, scoring_matrix, alignment_matrix)
        scoring_set[dummy] = score

    # compute un_normalized distribution
    scoring_list = list(scoring_set.values())
    scoring_dist = {}
    for dummy in scoring_list:
        if dummy != 0:
            count_score = scoring_list.count(dummy)
            scoring_dist[dummy] = count_score
    #compute normalized distribution
    # scoring_dist_norm = list(scoring_dist)
    sum_trials = sum(scoring_dist.values())
    x_data = list(scoring_dist.keys())
    y_data = []
    for dummy_i in list(scoring_dist.values()):
        y_data.append(dummy_i/float(sum_trials))
    plt.bar(x_data,y_data)
    plt.xlabel("max score local alignment")
    plt.ylabel("fraction of total trials")
    plt.title("Null Distribution for Hypothesis testing on 1000 Trials ")
    plt.show()
    print("scoring_dist = ", scoring_dist)
    return scoring_dist

num_trials = 1000
scoring_dist =  {45: 48, 49: 68, 44: 45, 48: 72, 52: 67, 46: 65, 50: 69, 40: 7, 47: 59, 51: 61, 43: 26, 59: 28, 60: 22,
                 42: 24, 55: 46, 41: 14, 54: 50, 61: 17, 65: 8, 53: 51, 58: 27, 66: 5, 57: 19, 63: 11, 62: 13, 73: 2, 56: 36,
                 67: 8, 78: 2, 80: 1, 77: 2, 64: 9, 71: 2, 39: 1, 68: 5, 69: 3, 90: 1, 91: 1, 70: 2, 76: 1, 75: 2}

def mean_standard_deviation_z_score(scoring_dist, num_trials):
    """

    :return:
    """
    #compute the mean
    scores_list = numpy.multiply(list(scoring_dist.keys()), list(scoring_dist.values()))
    print("sum(scoring_list) = ", sum(scores_list))

    mean = sum(scores_list)/float(num_trials)
    #compute the standard deviation
    sqrt_sum = 0
    for score in list(scoring_dist.keys()):
        num_score = scoring_dist.get(score)
        sqrt_sum += num_score*(score - mean)**2

    standard_dev = math.sqrt(sqrt_sum/float(num_trials))
    #compute z-score
    #score is the score of the local alignment for the human eyeless protein as computed is 875
    score_hum_ff = 875
    z_score = (score_hum_ff -mean)/standard_dev
    print("mean, standard_dev, z_score = ",mean, standard_dev, z_score)
    return mean, standard_dev, z_score

def edit_distance(seq_x, seq_y):
    """
    The edit distance corresponds to the minimum number of single character insertions, deletions, and substitutions
    that are needed to transform one string into another.
    The edit distance for two strings x and y can be expressed in terms of the lengths of the two strings and their corresponding similarity score
    as follows:∣x∣+∣y∣−score(x,y) where score(x,y) is the score returned by the global alignment of these two strings
    Determine the values for diag_score, off_diag_score, and dash_score such that the score from the resulting
    global alignment yields the edit distance when substituted into the formula above
    """
    # consider the case when seq_x = "a", seq_y = "a", edit_distance = 0 (no action needed),score(seq_x,seq_y) = 2 then diag_score = 2
    diag_score = 2
    # when seq_x = "a", seq_y = "", edit_distance = 1 (insert or delettion), score(seq_x,seq_y) = 0 then dash_score = 0
    dash_score = 0
    #when seq_x = "a", seq_y = "b", edit_distance = 1 (substitution), score(seq_x,seq_y) = 1 then off_diag_score = 1
    off_diag_score = 1
    return diag_score, off_diag_score, dash_score


def check_spelling(checked_word, dist, word_list):
    """
    iterates through word_list and returns the set of all words that are within edit distance dist of the string checked_word.
    :param checked_word:
    :param dist:
    :param word_list:
    :return:
    """
    result =[]
    for word in word_list:
        # print(word)
        diag_score, off_diag_score, dash_score = edit_distance(word, checked_word)
        alphabet = compute_alphabet(word, checked_word)
        scoring_matrix = build_scoring_matrix(alphabet, diag_score, off_diag_score, dash_score)
        alignment_matrix = compute_alignment_matrix(word, checked_word, scoring_matrix, global_flag=True)
        score_gl = compute_global_alignment(word, checked_word, scoring_matrix, alignment_matrix)
        edit_dist = len(word) + len(checked_word) - score_gl[0]
        if edit_dist <= dist:
            result.append(word)
    print(result)
    return result

word_list = read_word(CHECKT_LIST)
checked_word1 = "humble"
checked_word2 = "firefly"
dist1 = 1
dist2 = 2
# list1 = check_spelling(checked_word1, dist1, word_list)
# list2 = check_spelling(checked_word2, dist2, word_list)
