#**Assessment: Computing Alignments of Sequences and Applications to Genomics and Beyond**

#**Alignments of Sequences**

*Input:* Given an alphabet Σ and a scoring matrix M defined over Σ ∪ {′−′}, the dynamic programming method computed a score that measured the similarity of two sequences X and Y based on the values of this scoring matrix. In particular, this method involved computing an alignment matrix S between X and Y whose entry Sij​ scored the similarity of the substrings X[0…i−1] and Y[0…j−1].  
*Output:* Implement four functions. The first pair of functions will return matrices that we will use in computing the alignment of two sequences. The second pair of functions will return global and local alignments of two input sequences based on a provided alignment matrix. You will then use these functions in Application to analyze two problems involving comparison of similar sequences.  
#**Applications to Genomics**

**Comparing two proteins**  
For the genomics part of the Application, you will load several protein sequences and an appropriate scoring matrix. A substring of the eyeless protein of about 130 amino acids, known as the PAX domain, whose function is to bind specific sequences of DNA, is virtually identical between the human and fruit fly versions of eyeless.we compute the similarity between the human and fruit fly versions of the eyeless protein and see if we can identify the PAX domain.

- HumanEyelessProtein: http://storage.googleapis.com/codeskulptor-alg/alg_HumanEyelessProtein.tx
- FruitflyEyelessProtein: http://storage.googleapis.com/codeskulptor-alg/alg_FruitflyEyelessProtein.txt
- The scoring matrix PAM50: http://storage.googleapis.com/codeskulptor-alg/alg_PAM50.txt

Next, compute the local alignments of the sequences of HumanEyelessProtein and FruitflyEyelessProtein using the PAM50 scoring matrix.
consider the similarity of the two sequences in the local alignment computed in Question 1 to a third sequence. The file ConsensusPAXDomain contains a "consensus" sequence of the PAX domain; that is, the sequence of amino acids in the PAX domain in any organism.

- ConsensusPAXDomain: http://storage.googleapis.com/codeskulptor-alg/alg_ConsensusPAXDomain.txt

Load the file ConsensusPAXDomain. For each of the two sequences of the local alignment computed in Question 1, do the following:

    Delete any dashes '-' present in the sequence.
    Compute the global alignment of this dash-less sequence with the ConsensusPAXDomain sequence.
    Compare corresponding elements of these two globally-aligned sequences (local vs. consensus) and compute the percentage of elements in these two sequences that agree.  
To reiterate, you will compute the global alignments of local human vs. consensus PAX domain as well as local fruitfly vs. consensus PAX domain. Your answer should be two percentages: one for each global alignment. Enter each percentage below. Be sure to label each answer clearly and include three significant digits of precision.
s it likely that the level of similarity exhibited by the answers could have been due to chance? In particular, if you were comparing two random sequences of amino acids of length similar to that of HumanEyelessProtein and FruitflyEyelessProtein, would the level of agreement in these answers be likely? To help you in your analysis, there are 23 amino acids with symbols in the string ("ACBEDGFIHKMLNQPSRTWVYXZ". Include a short justification for your answer.

**Hypothesis testing**

Write a function generate_null_distribution(seq_x, seq_y, scoring_matrix, num_trials) that takes as input two sequences seq_x  and seq_y, a scoring matrix scoring_matrix, and a number of trials num_trials. This function should return a dictionary scoring_distribution that represents an un-normalized distribution generated by performing the following process num_trials times
    Generate a random permutation rand_y of the sequence seq_y using random.shuffle().
    Compute the maximum value score for the local alignment of seq_x and rand_y using the score matrix scoring_matrix.
    Increment the entry score in the dictionary scoring_distribution by one.

Use the function generate_null_distribution to create a distribution with 1000 trials using the protein sequences HumanEyelessProtein and FruitflyEyelessProtein (using the PAM50 scoring matrix).Next, create a bar plot of the normalized version of this distribution using plt.bar in matplotlib (or your favorite plotting tool).
To this end, we first compute the mean μ and the standard deviation σ of this distribution.

μ=(1/n)∑isi,

σ=sqrt((1/n)∑i(si−μ)**2)

where the values si​ are the scores returned by the n trials. If s is the score of the local alignment for the human eyeless protein and the fruitfly eyeless protein, the z-score zzz for this alignment is

z=(s−μ)/σ.

**Spelling correction**  

Given two strings, the edit distance corresponds to the minimum number of single character insertions, deletions, and substitutions that are needed to transform one string into another. In particular, if x and y are strings and a and b are characters, these edit operations have the form:

    Insert - Replace the string x+y by the string x+a+y.
    Delete - Replace the string x+a+y by the string x+y.
    Substitute - Replace the string x+a+y by the string x+b+y,

Not surprisingly, similarity between pairs of sequences and edit distances between pairs of strings are related. In particular, the edit distance for two strings xxx and yyy can be expressed in terms of the lengths of the two strings and their corresponding similarity score as follows:∣x∣+∣y∣−score(x,y) where score(x,y) is the score returned by the global alignment of these two strings using a very simple scoring matrix that can be computed using build_scoring_matrix.
Determine the values for diag_score, off_diag_score, and dash_score such that the score from the resulting global alignment yields the edit distance when substituted into the formula above. Be sure to indicate which values corresponds to which parameters. Finally, as a side note, be aware that there are alternative formulations of edit distance as a dynamic programming problem using different scoring matrices. For this problem, please restrict your consideration to the formulation used above.

http://storage.googleapis.com/codeskulptor-assets/assets_scrabble_words3.txt

To begin, load this list of 79339 words. Then, write a function check_spelling(checked_word, dist, word_list) that iterates through word_list and returns the set of all words that are within edit distance dist of the string checked_word.

Use your function check_spelling to compute the set of words within an edit distance of one from the string "humble" and the set of words within an edit distance of two from the string "firefly". (Note this is not "fruitfly".)

#**Outcomes:**  

**- Comparing two proteins:**  
Local alignment score: 875  
Humman = 'HSGVNQLGGVFVNGRPLPDSTRQKIVELAHSGARPCDISRILQVSNGCVSKILGRYYETGSIRPRAIGGSKPRVATPEVVSKIAQYKRECPSIFAWEIRDRLLSEGVCTNDNIPSVSSINRVLRNLASEK-QQ'  
Fruitfly = 'HSGVNQLGGVFVGGRPLPDSTRQKIVELAHSGARPCDISRILQVSNGCVSKILGRYYETGSIRPRAIGGSKPRVATAEVVSKISQYKRECPSIFAWEIRDRLLQENVCTNDNIPSVSSINRVLRNLAAQKEQQ'  
Length = 133  
Corresponding procentage similarities to  a third sequence ConsensusPAXDomain as follows:
* human_consensus_agree = 73% (97 characters to 133).
* fruitfly_consensus_agree = 70% ( 94 characters to 133).  
The probablity that 97 characters to match 133 elements sequence, whose elements choisen from a 23 character alphabet is less the 10^-100, so not likely the similarity between human vs. consensus PAX domain due to chance. 

**- Hypothesis testing:**  
From the generated grapic of the Null Hypothesis Statistic Test, which has the normal distribution,we can see that the vast majority of the distribution peak roughly centered at 50.   
From the formulas,the results given:

* mean = 52 
* standard_dev = 6.8
* z_score = 122  
We will assume that 99% of the scrores are within three standard deviations of the mean for this distribution. The acctual score of the Human/Fruityfly z_score is more than 100 standard deviations away from the mean of the distribution. If we assume that each multiple of three standard deviation reduces the likelt hood of this score arising randomly by 10^-2, the resulting proprability is around 10^-67
So the likelihood that the human sequence to math 73%  with any other random mutated sequences of fruityfly is almost zero.
 
**- Spelling correction:**  
Consider the case when seq_x = "a", seq_y = "a", edit_distance = 0 (no action needed),score(seq_x,seq_y) = 2 then diag_score = 2.  
when seq_x = "a", seq_y = "", edit_distance = 1 (insert or delettion), score(seq_x,seq_y) = 0 then dash_score = 0.  
when seq_x = "a", seq_y = "b", edit_distance = 1 (substitution), score(seq_x,seq_y) = 1 then off_diag_score = 1.  

list1 =['bumble', 'fumble', 'humble', 'humbled', 'humbler', 'humbles', 'humbly', 'jumble', 'mumble', 'rumble', 'tumble']
list2 = ['direly', 'finely', 'fireclay', 'firefly', 'firmly', 'firstly', 'fixedly', 'freely', 'liefly', 'refly', 'tiredly']

