# -*- coding: utf-8 -*-
"""
Created on Sat Aug  2 15:06:58 2025

@author: arnol
"""

import os
import csv
import pandas as pd
import numpy as np
import time


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
        A numerical representation of phonetic similarity between two phonemes. Differences in different features are weighed differently: place, voicing, manner contribute 1,2, and 4 to the difference, respectively. 
        
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



#%% 


#read file 
file = 'stops.csv'
df = pd.read_csv(file)

#go through the dataframe, and for each word and anadrome, compute 
#phonemic differences of first and last phoneme
for row in df.itertuples(): 
    #obtain word and anadrome transcription, stripping spaces
    word = row[5].strip()
    ana = row[6].strip()
    #get first and last phonemes of each 
    p1_word, p1_ana = word[1],ana[1]
    p3_word,p3_ana = word[-2],ana[-2]
    
    #compute phonemic difference of first and last phoneme 
    p1_diff = compute_phonemic_difference(p1_word,p1_ana)
    p3_diff = compute_phonemic_difference(p3_word,p3_ana)
    
    #compute total phonemic diff of word and anadrome
    word_diff = p1_diff + '0' + p3_diff
    #add the result to the dataframe
    index = row.Index #access index of that row
    df.loc[index, 'Phonemic Difference'] = word_diff
    
    
#%%

#now, we want to generate a list of controls with the same phonemic difference
#as the target word. first, create a list that will store controls 
controls_list = []

#create a simplified dictionary for CVC generation
filtered_phoneme_dict = {key: value for key, value in phoneme_dict.items() if key != "g"}

#iterate through dataframe again
for row in df.itertuples(): 
    #obtain word transcription
    transcription = row[5].strip()
    #obtain first, second, third phonemes 
    p1, p2, p3 = transcription[1],transcription[2:-2],transcription[-2]
    #given fixed vowel, iterate through the list of 
    exhaustive_cvcs =  ["/" + c1 + p2 + c2 + "/" for c1 in filtered_phoneme_dict for c2 in filtered_phoneme_dict]
    #find phonemic similarity for each control 
    cvc_similarity = [] #create a list to store this 
    for cvc in exhaustive_cvcs: 
        #get the first phoneme and 3rd phoneme of the cvc 
        control_p1, control_p3 = cvc[1],cvc[-2]
        p1_diff = compute_phonemic_difference(p1, control_p1)
        p3_diff = compute_phonemic_difference(p3, control_p3)
        word_diff = p1_diff + '0' + p3_diff 
        cvc_similarity.append(word_diff)
    controls_candidates = dict(zip(exhaustive_cvcs, cvc_similarity))
    #now we want to go through each entry of the controls dictionary
    #and keep only the ones that have the same phonemic difference from 
    #the main word as the anadrome does 
    phon_diff = row[14] #fetch word/ana phoneme difference
    controls_candidates_matching = {k: v for k, v in controls_candidates.items() if v == phon_diff}  
    #write the transcriptions themselves into the controls list 
    controls_list.append(controls_candidates_matching)
    # controls_list.append([key for key in controls_candidates_matching]) 
    
#visualize variable to determine controls 
# #%% 
# #write new data into csv 
# file = 'stops.csv'
# df.to_csv(file, index=False)

#%% Genenrate a graphical representations of words and their connections

import networkx as nx
import matplotlib.pyplot as plt


#Create a graph
G = nx.Graph()


#create a simplified dictionary for CVC generation
filtered_phoneme_dict = {key: value for key, value in phoneme_dict.items() if key != "g"}

#define vowel. nature of connections does not depend on vowel, only 
#whether each word is a node or not
vowel = 'ʌ'

#make exhaustive CVC list. this is our list of nodes  
nodes =  ["/" + c1 + vowel + c2 + "/" for c1 in filtered_phoneme_dict for c2 in filtered_phoneme_dict]
G.add_nodes_from(nodes)


#first loop --- iterate through all cvc trigrams    
for i, cvc in enumerate(nodes): 
    #obtain the first and third phoneme of the cvc 
    cvc_p1, cvc_p3 = cvc[1],cvc[-2]
    #first find anadrome, if it exists and is not itself 
    for cvc_ana in nodes: 
        #obtain first and third phoneme of cvc being used to compare
        cvc_ana_p1, cvc_ana_p3 = cvc_ana[1],cvc_ana[-2]
        #check whether anadrome condition satisfied 
        if '/' + cvc_ana_p3 + vowel + cvc_ana_p1 + '/' == cvc and cvc_ana != cvc: 
            G.add_edge(cvc, cvc_ana, connection='Anadrome')
            anadrome = cvc_ana
    #if an anadrome was not found, then the anadrome is the word itself
    if 'anadrome' not in locals(): 
        anadrome = cvc
    #compute phonemic difference between word and anadrome 
    #compute difference between first and last phonemes
    ana_p1,ana_p3 = anadrome[1],anadrome[-2]
    p1_diff = compute_phonemic_difference(cvc_p1,ana_p1)
    p3_diff = compute_phonemic_difference(cvc_p3,ana_p3)
    #compute total phonemic diff of word and anadrome
    ana_diff = p1_diff + '0' + p3_diff
    
    #now loop through remaining words looking for controls 
    remaining_cvc = nodes[i+1:]
    for cvc_ctrl in remaining_cvc: 
        #compute phonemic difference between word and 
        #control candidate
        control_p1, control_p3 = cvc_ctrl[1],cvc_ctrl[-2]
        p1_diff = compute_phonemic_difference(cvc_p1, control_p1)
        p3_diff = compute_phonemic_difference(cvc_p3, control_p3)
        word_diff = p1_diff + '0' + p3_diff 
        #check if equidistant from word as anadrome is 
        if word_diff == ana_diff and cvc_ctrl != anadrome: 
            #if yes, draw a control connection
            G.add_edge(cvc,cvc_ctrl, connection = 'Control')
    #once loop finishes, delete anadrome and continue 
    del anadrome

#%% Draw the graph

from mpl_toolkits.mplot3d import Axes3D
import random
import plotly.graph_objects as go

import plotly.io as pio

# Set renderer to open in the default web browser
pio.renderers.default = "browser"

# Define edge types and their colors
edge_colors = {
    'Anadrome': 'red',
    'Control': 'blue'
}

# Split the graph into connected components
subgraphs = [G.subgraph(c).copy() for c in nx.connected_components(G)]

# Loop through each subgraph and plot it
for idx, sg in enumerate(subgraphs):
    # Create random 3D positions for nodes
    pos = {node: (random.random(), random.random(), random.random()) for node in sg.nodes()}

    # Create node traces
    node_x = []
    node_y = []
    node_z = []
    node_text = []
    for node, coords in pos.items():
        x, y, z = coords
        node_x.append(x)
        node_y.append(y)
        node_z.append(z)
        node_text.append(str(node))

    node_trace = go.Scatter3d(
        x=node_x, y=node_y, z=node_z,
        mode='markers+text',
        marker=dict(size=6, color='lightgray'),
        text=node_text,
        textposition='top center',
        hoverinfo='text'
    )

    # Create edge traces for each type
    edge_traces = []
    for connection_type, color in edge_colors.items():
        edge_x = []
        edge_y = []
        edge_z = []
        for u, v, d in sg.edges(data=True):
            if d.get('connection') == connection_type:
                x0, y0, z0 = pos[u]
                x1, y1, z1 = pos[v]
                edge_x += [x0, x1, None]
                edge_y += [y0, y1, None]
                edge_z += [z0, z1, None]

        edge_trace = go.Scatter3d(
            x=edge_x, y=edge_y, z=edge_z,
            mode='lines',
            line=dict(color=color, width=4),
            hoverinfo='none',
            name=connection_type
        )
        edge_traces.append(edge_trace)

    # Combine all traces
    fig = go.Figure(data=[*edge_traces, node_trace])

    fig.update_layout(
        title=f"3D Subgraph {idx + 1}",
        showlegend=True,
        legend=dict(x=0, y=1),
        margin=dict(l=0, r=0, b=0, t=40),
        scene=dict(
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            zaxis=dict(visible=False)
        )
    )

    fig.show()
    
#%%
# # Define fricatives with attributes
# fricatives = set([
#     # Labiodental
#     Phoneme("f", "voiceless", "labiodental"),
#     Phoneme("v", "voiced", "labiodental"),

#     # Dental
#     Phoneme("θ", "voiceless", "dental"),  # theta
#     Phoneme("ð", "voiced", "dental"),    # eth

#     # Alveolar
#     Phoneme("s", "voiceless", "alveolar"),
#     Phoneme("z", "voiced", "alveolar"),

#     # Postalveolar
#     Phoneme("ʃ", "voiceless", "postalveolar"),  # "sh" sound
#     Phoneme("ʒ", "voiced", "postalveolar"),     # "zh" as in "measure"
#     ])

# fricative_dict = {p.symbol: p for p in fricatives}

# # Define the data as a dictionary
# data = {
#     'Word': ['/feɪs/', '/seɪv/'],
#     'Anadrome': ['/seɪf/', '/veɪs/'],
#     'Phonemic Difference': [0, 0]
# }

# df_fricatives = pd.DataFrame(data)

# #go through the dataframe, and for each word and anadrome, compute 
# #phonemic differences of first and last phoneme
# for row in df_fricatives.itertuples(): 
#     #obtain word and anadrome transcription, stripping spaces
#     word = row[1].strip()
#     ana = row[2].strip()
#     #get first and last phonemes of each 
#     p1_word, p1_ana = word[1],ana[1]
#     p3_word,p3_ana = word[-2],ana[-2]
    
#     #compute phonemic difference of first and last phoneme 
#     p1_diff = compute_phonemic_difference(p1_word,p1_ana)
#     p3_diff = compute_phonemic_difference(p3_word,p3_ana)
    
#     #compute total phonemic diff of word and anadrome
#     word_diff = p1_diff + '0' + p3_diff
#     #add the result to the dataframe
#     index = row.Index #access index of that row
#     df_fricatives.loc[index, 'Phonemic Difference'] = word_diff
            
    


    


