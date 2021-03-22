import yaml
import os
from pathlib import Path

currentDir = str(Path().absolute())  # app directory
possibleReplacementYml = "nucleotide_possibility.yml"
with open(possibleReplacementYml) as file:
    possible_repl_dict = yaml.load(file, Loader=yaml.FullLoader)

# print(possible_repl_dict)


def getProbePattern(probeName, seq_to_search):
    # Create pattern to search
    pattern = ""
    for idx, charVal in enumerate(seq_to_search):
        strAdd = "["

        # for G
        if charVal == "G":
            rep_list = possible_repl_dict[probeName][charVal]
            for el in rep_list:
                strAdd += el
        # for A
        elif charVal == "A":
            rep_list = possible_repl_dict[probeName][charVal]
            for el in rep_list:
                strAdd += el

        # for T
        elif charVal == "T":
            rep_list = possible_repl_dict[probeName][charVal]
            for el in rep_list:
                strAdd += el

        elif charVal == "C":
            rep_list = possible_repl_dict[probeName][charVal]
            for el in rep_list:
                strAdd += el

        strAdd += "]"
        pattern += strAdd

    return pattern

# def getProbePattern(probeName, seq_to_search):
#     # Create pattern to search
#     pattern = ""
#     if probeName == "Nucleocapsid_F":
#         for idx, charVal in enumerate(seq_to_search):
#             strAdd = "["
#             strAdd += charVal
#             # for A
#             if charVal == "A":
#                 rep_list = possible_repl_dict[charVal]
#                 for elidx, el in enumerate(rep_list):
#                     strAdd += el

#             # for T
#             elif charVal == "T":
#                 rep_list = possible_repl_dict[charVal]
#                 for elidx, el in enumerate(rep_list):
#                     strAdd += el

#             # for G
#             elif charVal == "G":
#                 rep_list = possible_repl_dict[charVal]
#                 for elidx, el in enumerate(rep_list):
#                     strAdd += el

#             # for C
#             elif charVal == "C":
#                 rep_list = possible_repl_dict[charVal]
#                 for elidx, el in enumerate(rep_list):
#                     strAdd += el

#             strAdd += "]"
#             pattern += strAdd

#     return pattern


def getProbePattern2(probeName, seq_to_search):
    # Create pattern to search
    pattern = ""
    if probeName == "Nucleocapsid_F":
        for idx, charVal in enumerate(seq_to_search):
            strAdd = "["
            strAdd += charVal
            # for A
            if charVal == "A":
                if idx in [17, 19]:
                    strAdd += "C"
                    strAdd += "M"

            # for T
            elif charVal == "T":
                strAdd += "C"
                strAdd += "A"

            # for G
            elif charVal == "G":
                if idx in [0, 1]:
                    strAdd += "A"
                    strAdd += "R"
                if idx in [2]:
                    strAdd += "C"
                if idx in [3]:
                    strAdd += "K"
                if idx in [18]:
                    strAdd += "T"

            # for C
            elif charVal == "C":
                if idx in [6, 15]:
                    strAdd += "Y"
                if idx in [6]:
                    strAdd += "T"

            strAdd += "]"
            pattern += strAdd

    elif probeName == "Nucleocapsid_R":
        for idx, charVal in enumerate(seq_to_search):
            strAdd = "["
            strAdd += charVal
            # for A
            if charVal == "A":
                if idx in [1, ]:
                    strAdd += "C"

            # for T
            elif charVal == "T":
                if idx in [21]:
                    strAdd += "A"

            # for G
            elif charVal == "G":
                if idx in [2, ]:
                    strAdd += "A"
                    strAdd += "T"
                    strAdd += "K"
                if idx in [10, 17, 22]:
                    strAdd += "T"
                    strAdd += "C"
                if idx in [18]:
                    strAdd += "K"

            # for C
            elif charVal == "C":
                if idx in [0]:
                    strAdd += "G"
                    strAdd += "Y"
                if idx in [3, ]:
                    strAdd += "S"
                if idx in [11, ]:
                    strAdd += "A"
                    strAdd += "T"
                if idx in [20, ]:
                    strAdd += "T"
                if idx in [11, ]:
                    strAdd += "Y"

            strAdd += "]"
            pattern += strAdd
    return pattern


def iterateText(foo):
    prevnl = -1
    while True:
        nextnl = foo.find('\n', prevnl + 1)
        if nextnl < 0:
            break
        yield foo[prevnl + 1:nextnl]
        prevnl = nextnl


def build_page(text):
    from xml.sax.saxutils import escape
    text = escape(text)
    return "<html><body>{0}</body></html>".format(text)
