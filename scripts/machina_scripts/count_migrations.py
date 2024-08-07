import sys
import numpy as np
from scipy.stats import entropy

def weird_division(n, d):
    return n / d if d else 0

def safe_entropy(x):
    if sum(x) == 0:
        h = np.nan
    else:
        h = entropy(x,base=3)
    return h

def tree_met_rate(metastasis_edges,non_metastasis_edges):
    total_edges = metastasis_edges + non_metastasis_edges
    rate = round(float(metastasis_edges / total_edges),5)
    return rate


### SET PARAMETERS ###
# Output transition counts or relative transition probabilities
transition_probs = False
ignore_self_migration = False
### END PARAMETERS

# read in big results file for all CP
# output one line per CP as well as one total line
# each line includes summary info as well as each cell in the transition table

col_dict = {}
lab_dict = {}
model = ""
migrations = []
cur_cp = ""

tissues = set()
cp = set(['all'])
with open(sys.argv[1],'r') as infile:
    lines = infile.readlines()
    for l in lines:
        l = l.strip()
        l = l.split(" ")
        cp.add(l[0])
        if l[1] == "color":
            tissues.add(l[2])


cp = sorted(list(cp))
tissues = sorted(list(tissues))

source_matrix_names = []
source_entropy_names = []

td = {}
for c in cp:
    td[c] = {}
    for t1 in tissues:
        td[c][t1] = {}
        if c == cp[0]:
            ent_name = "H_" + t1
            source_entropy_names.append(ent_name)
        for t2 in tissues:
            td[c][t1][t2] = 0
            # add transition names once
            if c == cp[0]:
                trans_name = t1 + ":" + t2
                source_matrix_names.append(trans_name)

met_edges = 0
non_met_edges = 0
all_met_edges = 0
all_non_met_edges = 0

# print header
out_header = ["CP,model,migrations,edges,TreeMetRate"] + source_entropy_names + source_matrix_names
print(*out_header,sep=",")

with open(sys.argv[1],'r') as infile:
    lines = infile.readlines()
    lines.append("END")
    for l in lines:
        l = l.strip()
        l = l.split(" ")
        cp = l[0]
        if cur_cp == "":
            cur_cp = cp
        #print(f'comparing {l} to last line {last}')
        if cp != cur_cp or l == ["END"]:
            # output CP summary
            # calculate entropy
            d = td[cur_cp]
            source_matrix = []
            source_count_matrix = []
            source_entropy = []
            all_tis_source_tot = 0
            all_source_probs = []
            for t1 in tissues:
                source_tot = 0
                target_tot = 0
                source_probs = []
                source_counts = []
                target_probs = []
                for t2 in tissues:
                    source_tot += d[t1][t2]
                    target_tot += d[t2][t1]
                    all_tis_source_tot += d[t1][t2]
                for t2 in tissues:
                     source_probs.append(weird_division(d[t1][t2],source_tot))
                     target_probs.append(weird_division(d[t2][t1],target_tot))
                     source_counts.append(d[t1][t2])
                hsource = safe_entropy(source_probs)
                htarget = safe_entropy(target_probs)
                source_matrix.append(source_probs)
                source_count_matrix.append(source_counts)
                source_entropy.append(hsource)
            for t1 in tissues:
                for t2 in tissues:
                    all_source_probs.append(weird_division(d[t1][t2],all_tis_source_tot))
            # Output probabilities
            if transition_probs:
                m = list(np.around(np.array(source_matrix),5))
            # Output counts instead
            else:
                m = list(np.around(np.array(source_count_matrix),5))

            source_entropy = list(np.around(np.array(source_entropy),5))
            source_entropy = ["NA" if np.isnan(x) else x for x in source_entropy]

            tmrate = tree_met_rate(met_edges,non_met_edges)
            edges = met_edges + non_met_edges
            outline = [cur_cp,model,len(migrations),edges,tmrate] + source_entropy + list(np.array(m).flatten())
            print(*outline,sep=",")
            #print(cur_cp,model,len(migrations),source_entropy[0],source_entropy[1],source_entropy[2],m[0][0],m[0][1],m[0][2],m[1][0],m[1][1],m[1][2],m[2][0],m[2][1],m[2][2],m,all_source_probs)
            col_dict = {}
            lab_dict = {}
            model = ""
            met_edges = 0
            non_met_edges = 0
            migrations = []
            cur_cp = cp
            if l == ['END']:
                break
        
        group = l[1]
        if group == "color":
            col_dict[l[3]] = l[2]
        elif group == "model":
            model = l[2]
        elif group == "migration":
            migrations.append([l[2],l[3]])
        elif group == "label":
            lab_dict[l[3]] = l[2]
        elif group == "tree":
            parent_col = lab_dict[l[2]]
            child_col = lab_dict[l[3]]
            parent_tis = col_dict[parent_col]
            child_tis = col_dict[child_col]
            if "XXX" not in l[3]:
                if parent_tis != child_tis:
                    met_edges += 1
                    all_met_edges += 1
                else:
                    non_met_edges += 1
                    all_non_met_edges += 1

            if "ASV" in l[3]:
                leaf = True
            else:
                leaf = False
            #if parent_tis == child_tis and leaf:
            #    pass
            if ignore_self_migration:
                if parent_tis != child_tis:
                    td[cp][parent_tis][child_tis] += 1
                    td['all'][parent_tis][child_tis] += 1
        
            else:
                td[cp][parent_tis][child_tis] += 1
                td['all'][parent_tis][child_tis] += 1

# print entropy for all
d = td['all']
source_matrix = []
source_count_matrix = []
source_entropy = []
all_tis_source_tot = 0
all_source_probs = []

for t1 in tissues:
    source_tot = 0
    target_tot = 0
    source_probs = []
    source_counts = []
    target_probs = []
    for t2 in tissues:
        source_tot += d[t1][t2]
        target_tot += d[t2][t1]
        all_tis_source_tot += d[t1][t2]
    for t2 in tissues:
         source_probs.append(weird_division(d[t1][t2],source_tot))
         target_probs.append(weird_division(d[t2][t1],target_tot))
         source_counts.append(d[t2][t1])
    hsource = safe_entropy(source_probs)
    htarget = safe_entropy(target_probs)
    source_matrix.append(source_probs)
    source_count_matrix.append(source_counts)
    source_entropy.append(hsource)
for t1 in tissues:
    for t2 in tissues:
        all_source_probs.append(weird_division(d[t1][t2],all_tis_source_tot))


s = list(np.around(np.array(source_entropy),5))
s = ["NA" if np.isnan(x) else x for x in s]
# use transition probabilities
if transition_probs:
    m = list(np.around(np.array(source_matrix),5))
# use transition counts
else:
    m = list(np.around(np.array(source_count_matrix),5))
tmrate = tree_met_rate(all_met_edges,all_non_met_edges)
all_edges = all_met_edges + all_non_met_edges
outline = ["all","NA","NA",all_edges,tmrate] + s + list(np.array(m).flatten())
# TO DO, order for 'all' columns being output is currently wrong
#print(*outline,sep=",")

