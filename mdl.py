from sys import argv, exit
from math import log
from copy import deepcopy
import numpy as np

pop = None
zipped_pop = None

class Group:
    def __init__(self,cal,*gene):
        self.gene_list = [g for g in gene]
        if cal:
            self.cal_bbwise_len()

    def make_copy(self):
        copy = Group(False,*self.gene_list)
        copy.gene_list = deepcopy(self.gene_list)
        copy.bb_d_model = self.bb_d_model
        copy.bb_d_data = self.bb_d_data
        return copy

    def append(self,grp):
        self.gene_list.extend(grp.gene_list)
        self.cal_bbwise_len()

    def cal_bbwise_len(self):
        combinations = pow(2,len(self.gene_list))
        self.bb_d_model = combinations - 1
        self.bb_d_data = 0
        bb_freq_list = [0 for _ in range(combinations)]
        
        selected_bits = deepcopy(zipped_pop[self.gene_list[0]])
        for i in range(1,len(self.gene_list)):
            selected_bits *= 2
            selected_bits += zipped_pop[self.gene_list[i]]

        for bits in selected_bits:
            bb_freq_list[bits] += 1

        # print(bb_freq_list)
        bb_freq_list = [item/len(pop) for item in bb_freq_list]
        # print(bb_freq_list)
        for item in bb_freq_list:
            self.bb_d_data += (item*log(item,2)) if item != 0 else 0
        # print(self.bb_d_data)


    def show(self):
        print(f'{len(self.gene_list)}',end=' ')
        for gene in self.gene_list:
            print(f'{gene}',end=' ')
        print("")


class MDL:
    def __init__(self):
        self.ell = len(pop[0])
        self.n = len(pop)
        self.group_list = [Group(True,i) for i in range(self.ell)]
        self.dscr_len = self.get_dscr_len(self.group_list)
    
    def get_dscr_len(self,group_list):
        d_model,d_data = 0,0
        for group in group_list:
            d_model += group.bb_d_model
            d_data += group.bb_d_data

        n = len(pop)
        d_model *= log(n,2)
        d_data *= -n

        return d_model + d_data

    def greedy_step(self):
        best_group_list = self.group_list
        best_dscr_len = self.dscr_len
        # print(f'base mdl = {self.dscr_len}')

        for m1 in range(len(self.group_list)):
            for m2 in range(m1+1,len(self.group_list)):
                new_group_list = [i.make_copy() for i in self.group_list]

                new_group_list[m1].append(new_group_list[m2])
                new_group_list.pop(m2)

                new_dscr_len = self.get_dscr_len(new_group_list)
                # print(new_dscr_len)
                if new_dscr_len < best_dscr_len:
                    # found better grouping
                    best_group_list = new_group_list
                    best_dscr_len = new_dscr_len

        # apply best grouping
        if self.group_list != best_group_list:
            self.group_list = best_group_list
            return True
        
        return False

    def greedy(self):
        while(self.greedy_step()):
            # print(len(self.group_list))
            pass

    def output(self):
        print(len(self.group_list))
        for group in self.group_list:
            group.show()


def read_file(path):
    l = []
    with open(path,'r') as f:
        while(True):
            line_str = f.readline()[:-1]
            if not len(line_str):
                break
            line = [int(bit) for bit in line_str]
            l.append(line)
    return l


if __name__ == '__main__':
    if len(argv)<2:
        print("usage: python3 mdl.py src_file")
        exit()

    pop = read_file(argv[1])
    zipped_pop = [list(i) for i in zip(*pop)]
    pop = np.array(pop)
    zipped_pop = np.array(zipped_pop)
    
    # print(zipped_pop[0])

    mdl = MDL()
    mdl.greedy()
    mdl.output()

