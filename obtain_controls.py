# -*- coding: utf-8 -*-
"""
Created on Sat Aug  2 15:06:58 2025

@author: arnol
"""



#define an abstract vowel representation. the consonants in controls do not depend 
#on the vowel. Feel free to replace with your favorite vowel. 
v = 'o'

#define phoneme object, which stores some information about the phoneme
class Phoneme:
    def __init__(self, symbol, manner, voicing, place):
        self.symbol = symbol  #e.g. "p"
        self.manner = manner #e.g. "fricative" or "plosive"
        self.voicing = voicing  #"voiced" or "voiceless"
        self.place = place  #e.g. "bilabial", "alveolar"

    def __repr__(self):
        return f"{self.symbol} ({self.voicing}, {self.place})"

#define a set of plosive phonemes with their attributes
plosives = set([
    Phoneme("p", "plosive", "voiceless", "bilabial"),
    Phoneme("b", "plosive", "voiced", "bilabial"),
    Phoneme("t", "plosive", "voiceless", "alveolar"),
    Phoneme("d", "plosive", "voiced", "alveolar"),
    Phoneme("k", "plosive", "voiceless", "velar"),
    Phoneme("ɡ", "plosive", "voiced", "velar"),  #correct IPA g (U+0261)
    Phoneme("g", "plosive", "voiced", "velar"),  #fallback if transcription uses ASCII (keyboard g)
])

#define a set of fricative phonemes with their attributes
fricatives = set([
    # Labiodental
    Phoneme("f", "fricative", "voiceless", "labiodental"),
    Phoneme("v", "fricative", "voiced", "labiodental"),

    # Dental
    Phoneme("θ", "fricative", "voiceless", "dental"),  # theta (unvoiced th)
    Phoneme("ð", "fricative", "voiced", "dental"),    # eth (voiced th)

    # Alveolar
    Phoneme("s", "fricative", "voiceless", "alveolar"),
    Phoneme("z", "fricative", "voiced", "alveolar"),

    # Postalveolar
    Phoneme("ʃ", "fricative", "voiceless", "postalveolar"),  # "sh" sound
    Phoneme("ʒ", "fricative", "voiced", "postalveolar"),     # "zh" as in "measure"
    ])

f_and_p = plosives | fricatives

#create dictionaries to access phoneme information from a string representation
plosive_dict = {p.symbol: p for p in plosives}
fricative_dict = {f.symbol: f for f in fricatives}
phoneme_dict = {p.symbol: p for p in f_and_p}

#define a function which takes letters as args, and returns a measure
#of how different they are.
#1 = diff voicing; 2 = diff place; 4 = diff manner (plosive vs fricative)
def compute_phonemic_difference(p1: str, p2: str) -> str: 
    """
    Calculates the difference beteween two phonemes in phonetic space. Inputs must be Unicode characters. 
    
    Parameters
    ----------
    p1 : str
        The first phoneme that you'd like to compare.
    p2 : str
        The second phoneme that you'd like to compare.

    Returns
    -------
    str
        A numerical representation of phonetic similarity between two phonemes. Differences in different features are weighed differently: place, voicing, manner contribute 1, 2, and 4 to the difference, respectively. 
        
    Examples:
    -------
    >>>compute_phonemic_difference('b','s')
    '7'
    
    >>>compute_phonemic_difference('t','d')
    '2'
    """
    
    diff = 0
    #check voicing difference
    if phoneme_dict.get(p1).voicing != phoneme_dict.get(p2).voicing: 
        diff +=2
    #check place difference
    if phoneme_dict.get(p1).place != phoneme_dict.get(p2).place: 
        diff += 1
    #check manner difference
    if phoneme_dict.get(p1).manner != phoneme_dict.get(p2).manner: 
        diff +=4 
    diff = str(diff)
    return diff


def compute_word_difference(word1: str, word2: str) -> str:
    """
    Calculates the difference between two tokens in phonetic space, specifically checking the starting and ending consonants. Inputs must be CVCs in IPA format, where both consonants exist in the phoneme set, and have the same vowel, e.g. /bʌt/ and /tʌb/.

    Parameters
    ----------
    word1 : str
        The first token of two that you want to compare.
    word2 : str
        The second token of two that you want to compare.

    Returns
    -------
    str
         A numerical representation of similarity in phonetic space between two tokens. Uses the compute_phonemic_similarity function.

    Examples
    --------
    >>>compute_word_difference('/bos/', '/sob/')
    '707'
    
    >>>compute_word_difference('/tod/', '/dot/') 
    '202'
    
    >>>compute_word_difference('/bos/', '/dot/')
    '104'
    """
    
    diff1 = compute_phonemic_difference(word1[1],word2[1])
    diff2 = compute_phonemic_difference(word1[-2],word2[-2])
    return diff1 + '0' + diff2

#Now, our goal is to create a dictionary where we can look up a CVC word and immediately access 
#appropriate control words. 

#generate an exhaustive list of SVS and SVF structures. 
#create simplified dictionaries to avoid redundancy
filtered_phoneme_dict = {key: value for key, value in phoneme_dict.items() if key != "g"}
filtered_plosive_dict = {key: value for key, value in plosive_dict.items() if key != "g"}


#create an exhaustive list of CVCs (fixing a vowel) 
all_svs =  ["/" + c1 + v + c2 + "/" for c1 in filtered_plosive_dict for c2 in filtered_plosive_dict]
all_svf = ["/" + c + v + f + "/" for c in filtered_plosive_dict for f in fricative_dict]
all_fvs = ["/" + f + v + c + "/" for c in filtered_plosive_dict for f in fricative_dict]
all_cvc = all_svs + all_svf + all_fvs

controls = []
#loop through the list of cvcs
for index, item in enumerate(all_cvc): 
    anadrome = item[::-1] #find the anadrome for the item 
    diff = compute_word_difference(item, anadrome) #compute phonemic difference between word/anadrome
    control_candidates = [] #create a list in which to store controls for a single token
    #loop through cvcs a second time to find control candidates
    for cvc in all_cvc: 
        #compute difference between word and cvc 
        control_diff = compute_word_difference(item, cvc)
        #if cvc is equidistant from word as anadrome is: 
        if control_diff == diff: 
            #the word becomes a control
            control_candidates.append((cvc, control_diff))
    #add to master control list 
    controls.append(control_candidates)

#create a dictionary from which to access controls for each word
word_control_dict = dict(zip(all_cvc, controls))






    





