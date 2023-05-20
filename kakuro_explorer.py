import matplotlib.pyplot as plt
import matplotlib.collections as collections
import numpy as np
import sys


b = 0.025
FONTSIZE_VAL = 10
FONTSIZE_HINT = 8

def fill_one_empty(ax, row, col):

    row = row - 1
    col = col - 1
    
    s_c = [col, col+1, col+1, col]
    s_r = [row, row, row+1, row+1]
    
    ax.fill(s_c, s_r, 'k')

    return ax

def fill_one_empty_g(ax, row, col):

    row = row - 1
    col = col - 1
    
    s_c = [col, col+1, col+1, col]
    s_r = [row, row, row+1, row+1]
    
    ax.fill(s_c, s_r, 'g')

    return ax

def fill_one_hint(ax, row, col, hint1, hint2):

    s_hint1 = str(hint1)
    s_hint2 = str(hint2)
    
    ax.fill([col-1+b,col,col], [row-1,row-1,row-b], 'k') 
    ax.fill([col-1,col-1,col-b], [row-1+b,row,row], 'k') 

    
    ax.text(col-1+0.5, row+-1+0.35, s_hint1, color='w', fontsize=FONTSIZE_VAL, fontweight='bold')
    ax.text(col-1+0.2, row+-1+0.8, s_hint2, color='w', fontsize=FONTSIZE_VAL, fontweight='bold')
    
    return ax

def fill_hints(ax, row, col, hintS):

    hint = hintS[0]
    _type = hintS[1]

    ax.fill([col-1+b,col,col], [row-1,row-1,row-b], 'tab:grey')
    ax.fill([col-1,col-1,col-b], [row-1+b,row,row], 'tab:grey')

    if(_type == "down"):
        ax.text(col-1+0.1, row+-1+0.9, hint[0], color='w', fontsize=FONTSIZE_HINT, fontweight='bold')
    if(_type == "right"):
        ax.text(col-1+0.6, row+-1+0.25, hint[0], color='w', fontsize=FONTSIZE_HINT, fontweight='bold')
    if(_type == "combined"):
        ax.text(col-1+0.1, row+-1+0.9, hint[0], color='w', fontsize=FONTSIZE_HINT, fontweight='bold')
        ax.text(col-1+0.6, row+-1+0.25, hint[1], color='w', fontsize=FONTSIZE_HINT, fontweight='bold')
        
def write_value(ax, row, col, val):

    row = row - 1
    col = col -1

    ax.text(col + 0.5, row + 0.5, val, color='k', fontsize=FONTSIZE_VAL, fontweight='bold', ha='center', va='center')
    
    
    return ax

def write_mark(ax, row, col):

    row = row - 1
    col = col -1

    ax.text(col + 0.5, row + 0.5, "?", color='g', fontsize=FONTSIZE_VAL, fontweight='bold', ha='center', va='center')


def write_value_s(ax, row, col, val):

    row = row - 1
    col = col -1

    ax.text(col + 0.5, row + 0.5, val, color='k', fontsize=FONTSIZE_VAL, fontweight='bold', ha='center', va='center')
    

    return ax

def create_matrix(m, n):

    print("m:", m, "n:",n)
    
    fig, ax = plt.subplots()
    ax.invert_yaxis()
    

    for i in range(m+1):
        ax.plot([0,m], [i,i], 'k-')
        ax.plot([i,i], [0,m], 'k-')

    for j in range(m, n+1):
        
        ax.plot([m,n], [j-m,j-m], 'k-')
        ax.plot([j,j], [0,m], 'k-')
    

    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
        
    return fig, ax

def assess_hint(hint):

    hint = hint / 10
    
    if(hint % 100 == 0):
        return [[int(hint/100)], "right"]
        
    else:
        div = int(hint / 100)
        rem = int(hint % 100)

        if(div == 0 and rem != 0):
            return [[rem], "down"]
        if(div != 0 and rem != 0):
            return [[rem, div], "combined"] ##[0] down [1] right


def read_file(_file):

    reader = open(_file, "r")
    readlines = reader.readlines()

    m, n = readlines[0].strip("\n").split(" ")
    
    m = int(m)
    n = int(n)
    
    board = [[] for i in range(m)]
    print("board: ", board)
    for i in range(1, len(readlines)):
        vals = readlines[i].strip("\n").split(" ")
        vals = [x for x in vals if x != ""]
        print("i: ", i , "vals: ", vals)
        for val in vals:
            board[i-1].append(int(val))

    return board, m, n
        
class Kakuro:
    
    def __init__(self, _file):

        global FONTSIZE
        
        self. board, self.m, self.n = read_file(_file)
        self.fig, self.ax = create_matrix(self.m , self.n)
        if(self.m > 20 or self.n > 20):
            FONTSIZE_VAL = 7
            FONTSIZE_HINT = 5
        
    def print_matrix(self):
        plt.show()


    def generate_board(self):
        
        boardSol = self.board
        
        i = 0
        for row in boardSol:
            j = 0
            for item in row:
                
                if(item == -1):
                    fill_one_empty(self.ax, i+1, j+1)
                elif(item > 0 and item <= 9):
                    write_value(self.ax, i+1, j+1, item)
                elif(item == -2):
                    write_mark(self.ax, i+1, j+1)
                else:
                    fill_hints(self.ax, i+1, j+1, assess_hint(item))
                j = j + 1
            i = i + 1
        plt.show()
    
if __name__ == "__main__":

    _file = sys.argv[1]
    k = Kakuro(_file)
    k.generate_board()
